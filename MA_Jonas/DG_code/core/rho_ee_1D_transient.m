clc, clearvars
%% Matrizen und Parameter
N = 1; % Polynomgrad des DG-Ansatzes, N+1 Knoten
K_werte = [12, 24, 48]; %, 96]; % Unterteilung der X-Richtung in 12 DG-Elemente
N_kx_werte = [15, 31, 61]; %, 121]; % Unterteilung der xi-Richtung in 15 FV-
                                % Zellen, ungerade für Mittelwert bei 0
n_Fehler = zeros(size(K_werte));
j_Fehler = zeros(size(K_werte));
N_Fehler = zeros(size(K_werte));
herm_Fehler = zeros(size(K_werte));

% Schleife für Gitterkonvergenz
for q = 1:length(K_werte)
    K = K_werte(q);
    N_kx = N_kx_werte(q);
    L_X = 4e-6; % Teststrecke von 4 µm
    d_X = L_X/K; % Breite eines DG-Elements
    d_kx = pi/L_X; % Breite eines DG-Elements
    
    m = (N_kx-1)/2;
    
    kx_werte = (-m:m)*d_kx;


    r = JacobiGL(0,0,N); % Berechnung der DG-Knoten: Vektor r=(-1,1)^T
    V = Vandermonde1D(N,r); % Vandermondematrix
    D = Dmatrix1D(N,r,V); % Ableitungsmatrix
    M = inv(V*V'); % Massematrix
    S = M*D; % Steifigkeitsmatrix
    M_k = (d_X/2)*M; % Massematrix mit physikalischer Elementbreite (Für 
                    % Steifigkeitsmatrix kürzt sich alles und es kann 
                    % gewöhnlich S verwendet werden.)
    t_End = 0.6e-12; % Festlegung der Endzeit (0.4ps für 0.4µm), für RK4 später
    
    
    % Berechnung der X-Knoten im jeweiligen DG-Element und der xi-Zellenmittelpunkte
    X = zeros(2,K);
    for k=1:K
        X(:,k) = 0.5*d_X*r + 0.5*(k*d_X + (k-1)*d_X);
    end
    
    %% Dirac-Transportmatrix
    vF = 1e6;
    A = vF*[0 1; 1 0];
    I = eye(2);
    A_X = 0.5*(kron(I,A) + kron(A.',I));
    A_xi = kron(I,A) - kron(A.',I);
    
    % Eigenwertzerlegung für Upwind
    % Behandlung X
    [EV_X,EW_X] = eig(A_X);
    lambda_X = diag(EW_X);
    A_X_plus  = EV_X*diag(max(lambda_X,0))*EV_X';
    A_X_minus = EV_X*diag(min(lambda_X,0))*EV_X';

    
    
    %% Verifikation des freien Transports eines e-Wellenpakets in Graphen.
    % Dazu als Wellenfunktion ein Gaußpaket (oder ebene Welle)
    x0 = 1.2e-6; % Position des Paketzentrums [m]
    sigma = 0.4e-6; % Paketbreite [m]
    g = @(x) exp(-0.5*(mod(x-x0+L_X/2,L_X)-L_X/2).^2/sigma^2); % Wellenpaket mit 
                    % mod zur Berücksichtigung periodischer Randbedingungen. 
    chi_0 = [1;0];
    rho_anf = chi_0*chi_0';
    
    U = zeros(4,2,K,N_kx); % Zusammenbauen des Vektors u
    for k = 1:K
        for j = 1:N_kx
            for n = 1:2
                rho_xi = @(xi) g(X(n,k)+xi/2).*g(X(n,k)-xi/2);
                f = 1/(2*pi)*integral(@(xi) rho_xi(xi) ...
                    .*exp(-1i*kx_werte(j)*xi), -L_X,L_X)*rho_anf;
                U(:,n,k,j) = f(:);
            end
        end
    end
    U_anfang = U; % für Vektorgestalt
    
    
    %% Zeitentwicklung: explizites Runge-Kutta-Verfahren 4. Ordnung
    omega_max = 2*vF*max(abs(kx_werte)); % Oszillation der Komponenten der 
                        % Wignerfunktion
    dt_max = 0.05*min(d_X/vF,1/omega_max); % Schnellste Oszillation oder 
                        % Zeit, die ein Teilchen zum Durchqueren eines
                        % DG-Elements benötigt
    n_steps = ceil(t_End/dt_max); % Anzahl der Zeitschritte insgesamt
    dt = t_End/n_steps; % Schrittweite (h_k im Skript)
    % Definition des zu lösenden Problems und Lösungsschleife
    rechts = @(W) DG_FV_rechts(W, S, M_k, kx_werte, A_X, A_X_plus, ...
        A_X_minus,A_xi);
    for it = 1:n_steps
        k1 = rechts(U);
        k2 = rechts(U + 0.5*dt*k1);
        k3 = rechts(U + 0.5*dt*k2);
        k4 = rechts(U + dt*k3);
        U = U + (dt/6)*(k1+2*k2+2*k3+k4);
    end
end



%% Ausgabe
xPlot = X(:) * 1e6;

rho11 = reshape(real(squeeze(U(1,:,:,:))), [], numel(kx_werte));
rho21 = reshape(real(squeeze(U(2,:,:,:))), [], numel(kx_werte));
rho12 = reshape(real(squeeze(U(3,:,:,:))), [], numel(kx_werte));
rho22 = reshape(real(squeeze(U(4,:,:,:))), [], numel(kx_werte));

figure

subplot(2,2,1)
mesh(kx_werte, xPlot, rho11)
xlabel('k_x [1/m]')
ylabel('X [\mum]')
zlabel('Re(\rho_{11})')
title('\rho_{11}')
grid on

subplot(2,2,2)
mesh(kx_werte, xPlot, rho12)
xlabel('k_x [1/m]')
ylabel('X [\mum]')
zlabel('Re(\rho_{12})')
title('\rho_{12}')
grid on

subplot(2,2,3)
mesh(kx_werte, xPlot, rho21)
xlabel('k_x [1/m]')
ylabel('X [\mum]')
zlabel('Re(\rho_{21})')
title('\rho_{21}')
grid on

subplot(2,2,4)
mesh(kx_werte, xPlot, rho22)
xlabel('k_x [1/m]')
ylabel('X [\mum]')
zlabel('Re(\rho_{22})')
title('\rho_{22}')
grid on

%% Dichtematrix zum Zeitpunkt t = 0
rho11_0 = reshape(real(squeeze(U_anfang(1,:,:,:))), [], numel(kx_werte));
rho21_0 = reshape(real(squeeze(U_anfang(2,:,:,:))), [], numel(kx_werte));
rho12_0 = reshape(real(squeeze(U_anfang(3,:,:,:))), [], numel(kx_werte));
rho22_0 = reshape(real(squeeze(U_anfang(4,:,:,:))), [], numel(kx_werte));

figure

subplot(2,2,1)
mesh(kx_werte, xPlot, rho11_0)
xlabel('k_x [1/m]')
ylabel('X [\mum]')
zlabel('Re(\rho_{11})')
title('\rho_{11}')
grid on

subplot(2,2,2)
mesh(kx_werte, xPlot, rho12_0)
xlabel('k_x [1/m]')
ylabel('X [\mum]')
zlabel('Re(\rho_{12})')
title('\rho_{12}')
grid on

subplot(2,2,3)
mesh(kx_werte, xPlot, rho21_0)
xlabel('k_x [1/m]')
ylabel('X [\mum]')
zlabel('Re(\rho_{21})')
title('\rho_{21}')
grid on

subplot(2,2,4)
mesh(kx_werte, xPlot, rho22_0)
xlabel('k_x [1/m]')
ylabel('X [\mum]')
zlabel('Re(\rho_{22})')
title('\rho_{22}')
grid on

sgtitle('Wigner-Matrix zum Zeitpunkt t = 0')

%% Lokale Funktion: rechte Seite fuer alle DG-Elemente und FV-Zellen
function dUdt = DG_FV_rechts(U, S, M_k, kx_werte, A_X, A_X_plus, A_X_minus, ...
    A_xi)
    K = size(U,3);
    N_kx = size(U,4);
    dUdt = zeros(size(U));
    % Schleife über alle DG-Elemente
    for k = 1:K
        k_links = mod(k-2,K)+1; % periodischer linker Nachbar
        k_rechts = mod(k,K)+1; % periodischer rechter Nachbar
        % Schleife über FV-Zellen
        for j = 1:N_kx
            U_mitte = U(:,:,k,j);  % 4x2-Matrix, die die Dichtematrixwerte 
                                   % eines DG-Elements zu eier FV-Zelle
                                   % enthält

            % DG-Randfluesse in X
            F_links = A_X_plus*U(:,2,k_links,j) + A_X_minus*U_mitte(:,1);
            F_rechts = A_X_plus*U_mitte(:,2) + A_X_minus*U(:,1,k_rechts,j);

            % Schwache DG-Formulierung
            F = A_X*U_mitte;
            linke_Seite = (F*S - F_rechts*[0 1] + F_links*[1 0])/M_k;

            % FV-Flussdifferenz in xi
            G = 1i*kx_werte(j)*A_xi*U_mitte;

            % Zeitliche Ableitung: dU/dt = DG - G
            dUdt(:,:,k,j) = linke_Seite-G;
        end
    end
end