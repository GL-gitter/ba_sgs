clc, clearvars
%% Matrizen und Parameter
N = 1; % Polynomgrad des DG-Ansatzes, N+1 Knoten
K_werte = [12, 24, 48]; %, 96]; % Unterteilung der X-Richtung in 12 DG-Elemente
N_xi_werte = [15, 31, 61]; %, 121]; % Unterteilung der xi-Richtung in 15 FV-
                                % Zellen, ungerade für Mittelwert bei 0
n_Fehler = zeros(size(K_werte));
j_Fehler = zeros(size(K_werte));
N_Fehler = zeros(size(K_werte));
herm_Fehler = zeros(size(K_werte));

hbar = 1.054571817e-34; % [Js]
e = 1.602176634e-19; % [C]
U0 = 2e-3*e; % Barrierenhöhe: 2 meV [J]
x_B = 2.2e-6; % Barrierenposition [m]
w_B = 0.25e-6; % Barrierenbreite [m]

U_pot = @(x) U0*exp(-0.5*((x-x_B)/w_B).^2);

% D_F = 1e-3*e/1e-6; % Potentialanstieg [J/m]
% U_pot = @(x) D_F*x;

% Schleife für Gitterkonvergenz
for q = 1:length(K_werte)
    K = K_werte(q);
    N_xi = N_xi_werte(q);
    L_X = 4e-6; % Teststrecke von 4 µm
    d_X = L_X/K; % Breite eines DG-Elements
    d_xi = 2*L_X/N_xi; % da xi=x_1-x_2 und x_1,x_2 element [0,L] sind, ist xi 
                       % element [-L,L]. Mit der Anzahl an Zellen ergibt sich 
                       % 2L/N
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
    xi = zeros(1,N_xi);
    for j=1:N_xi
        xi(j) = -L_X + (j - 0.5)*d_xi;
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
    % Behandlung xi
    [EV_xi,EW_xi] = eig(A_xi);
    lambda_xi = diag(EW_xi);
    A_xi_plus  = EV_xi*diag(max(lambda_xi,0))*EV_xi';
    A_xi_minus = EV_xi*diag(min(lambda_xi,0))*EV_xi';
    
    
    %% Implementation des von S_pot
    S_pot = zeros(2,K,N_xi);
    for k = 1:K
        for j = 1:N_xi
            x_1 = X(:,k)+xi(j)/2;
            x_2 = X(:,k)-xi(j)/2;

            S_pot(:,k,j) = (U_pot(x_1) - U_pot(x_2))/(1i*hbar);
        end
    end


    %% Verifikation des freien Transports eines e-Wellenpakets in Graphen.
    % Dazu als Wellenfunktion ein Gaußpaket (oder ebene Welle)
    x0 = 1.2e-6; % Position des Paketzentrums [m]
    sigma = 0.4e-6; % Paketbreite [m]
    g = @(x) exp(-0.5*((x-x0)/sigma).^2); % Wellenpaket mit 
                    % mod zur Berücksichtigung periodischer Randbedingungen. 
    chi_0 = [1;1]/sqrt(2);
    rho_anf = chi_0*chi_0';
    
    U = zeros(4,2,K,N_xi); % Zusammenbauen des Vektors u
    for k = 1:K
        for j = 1:N_xi
            for n = 1:2
                x1 = X(n,k) + xi(j)/2;
                x2 = X(n,k) - xi(j)/2;
                rho = g(x1)*g(x2)*rho_anf;
                U(:,n,k,j) = rho(:);
            end
        end
    end
    U_anfang = U; % für Vektorgestalt


    %% xi-Abhängigkeit von rho zum Zeitpunkt t=0
    % Für die Abhängigkeit wird ein fixer X-Wert benötigt. Dazu wird ein Punkt
    % nahe des Mittelpunkts gewählt. 
    % DG-Knoten nahe dem Zentrum des Gaußpakets auswählen:
    [~, idx] = min(abs(X(:)-x0)); % Index nahe x0 kriegen, Abstand irrelevant
    [n0,k0] = ind2sub(size(X),idx); % n0 und k0 zu dem Index speichern
    % Spur der Dichtematrix bei festem X für alle xi-Werte
    rho_xi = squeeze(U_anfang(1,n0,k0,:) + U_anfang(4,n0,k0,:) );
    rho_xi_ana = g(X(n0,k0)+xi/2) .* g(X(n0,k0)-xi/2);
    % Plot
    figure;
    plot(xi*1e6, real(rho_xi), 'bo', 'LineWidth', 1.5, 'DisplayName', 'Numerisch');
    hold on;
    plot(xi*1e6, real(rho_xi_ana), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Analytisch');
    xlabel('\xi [\mum]');
    ylabel('\rho(X,\xi,0)');
    title('Anfängliche xi-Abhängigkeit bei X nahe x_0');
    legend('Location','best');
    grid on;
    
    
    %% Zeitentwicklung: explizites Runge-Kutta-Verfahren 4. Ordnung
    dt_max = 0.05*min(d_X/vF,d_xi/(2*vF)); % CFL Bedingung maximal möglicher 
                % Zeitschritt, c=0.05 (Testwert). Minimum beider Werte, damit 
                % keine Zelle ggf. übersprungen wird.
    n_steps = ceil(t_End/dt_max); % Anzahl der Zeitschritte insgesamt
    dt = t_End/n_steps; % Schrittweite (h_k im Skript)
    % Definition des zu lösenden Problems und Lösungsschleife
    rechts = @(W) DG_FV_rechts(W, S, M_k, d_xi, A_X, A_X_plus, A_X_minus, ...
        A_xi_plus, A_xi_minus, S_pot);
    for it = 1:n_steps
        k1 = rechts(U);
        k2 = rechts(U + 0.5*dt*k1);
        k3 = rechts(U + 0.5*dt*k2);
        k4 = rechts(U + dt*k3);
        U = U + (dt/6)*(k1+2*k2+2*k3+k4);
    end


    % %% Vergleich mit analytischer Dichtematrix
    % rho_Fehler_summe = 0;
    % rho_Norm_summe = 0;
    % for k = 1:K
    %     for j = 1:N_xi
    %         for n = 1:2
    %             x1 = X(n,k) + xi(j)/2;
    %             x2 = X(n,k) - xi(j)/2;
    %             % Analytische freie Dichtematrix
    %             rho_frei = g(x1-vF*t_End)*g(x2-vF*t_End)*rho_anf;
    %             % Phasenfaktor des linearen Potentials
    %             phase = exp(-1i*D_F*xi(j)*t_End/hbar);
    %             % Analytische Lösung mit Potential
    %             rho_ana = phase*rho_frei;
    %             % Numerische Lösung
    %             rho_num = reshape(U(:,n,k,j),2,2);
    %             % Fehler aufsummieren
    %             rho_Fehler_summe = rho_Fehler_summe ...
    %                 + norm(rho_num-rho_ana,'fro')^2;
    %             rho_Norm_summe = rho_Norm_summe ...
    %                 + norm(rho_ana,'fro')^2;
    %         end
    %     end
    % end
    % rho_relFehler = sqrt(rho_Fehler_summe/rho_Norm_summe);
    % fprintf('Relativer Dichtematrixfehler: %.3e\n',rho_relFehler);
    
    
    %% Überprüfung: Teilchendichte, globale Teilchenzahlerhaltung, Hermitizität
    % und Vergleich mit der analytischen Lösung
    % 1) Teilchendichte
    j0 = (N_xi+1)/2; % Index der FV-Zelle, dessen Mittelpunkt bei xi=0 liegt
    n_anfang = squeeze(U_anfang(1,:,:,j0) + U_anfang(4,:,:,j0)); % n zu t=0
    n_ende = squeeze(U(1,:,:,j0) + U(4,:,:,j0)); % n im Endzustand
    % 2) Globale Teilchenzahlerhaltung
    N_anfang = (d_X/2)*sum(n_anfang(:)); % Trapezregel zur Auswertung des Integrals
    N_ende = (d_X/2)*sum(n_ende(:));
    N_relFehler = abs(N_ende-N_anfang)/N_anfang;
    % 3) Hermitizitaet: rho(X,xi) = rho(X,-xi)^dagger
    hermFehler = 0;
    for j = 1:N_xi
        for k = 1:K
            for n = 1:2
                % 2x2 Matrixdarstellung mit reshape
                rho_xi_plus = reshape(U(:,n,k,j),2,2); 
                rho_xi_minus = reshape(U(:,n,k,N_xi+1-j),2,2); % Zugreifen auf
                        % den Spiegelpunkt von xi_j
                % Fehler mit Frobenius Norm
                hermFehler = max(hermFehler,norm(rho_xi_plus-rho_xi_minus','fro'));
            end
        end
    end
    % 4) Analytischer Vergleich
    % Für die Teilchendichte
    n_ana = g(X-vF*t_End).^2;
    n_relFehler = norm(n_ende(:)-n_ana(:))/norm(n_ana(:));
    % Für den Strom
    j = squeeze(vF*(U(2,:,:,j0)+U(3,:,:,j0)));
    j_ana = vF*(g(X-vF*t_End).^2);
    j_relFehler = norm(j(:)-j_ana(:))/norm(j_ana(:));

    % Fehler
    n_Fehler(q) = n_relFehler;
    j_Fehler(q) = j_relFehler;
    N_Fehler(q) = N_relFehler;
    herm_Fehler(q) = hermFehler;
    
    % Überprüfung von Klein-Tunneln
    n_plus = real(1/2*(n_ende+j/vF));
    n_minus = real(1/2*(n_ende-j/vF));
    % Relativer Beitrag des linkslaufenden Anteils
    rel_L_Anteil = sum(abs(n_minus(:)))/sum(abs(n_ende(:)));
    fprintf('Relativer linkslaufender Anteil: %.3e\n', rel_L_Anteil);
    % Plot
    figure;
    plot(X(:)*1e6, n_plus(:), 'b.-', ...
        'DisplayName','Rechtslaufend'); hold on;

    plot(X(:)*1e6, n_minus(:), 'r.-', ...
        'DisplayName','Linkslaufend');

    xlabel('X [\mum]');
    ylabel('Dichte');
    title('Rechts- und linkslaufende Anteile');
    legend('Location','best');
    grid on;


    fprintf('K = %d, N_xi = %d, Dichtefehler = %.3e\n', ...
        K, N_xi, n_relFehler);
end

% Für die verschiedenen K verändert sich auch die Gitterweite: h=d_X=L_X/K
h = L_X./K_werte;
% Beobachtete Konvergenzordnung. Es ist epsilon1=Ch_1^p und
% epsilon2=ch_2^p, also p=ln(epsilon1/epsilon2)/ln(h_1/h_2)
% Teilchendichte
p_n = log(n_Fehler(1:end-1)./n_Fehler(2:end)) ./ log(h(1:end-1)./h(2:end));
% Strom
p_j = log(j_Fehler(1:end-1)./j_Fehler(2:end)) ./ log(h(1:end-1)./h(2:end));


%% Ausgabe
fprintf('Relativer Teilchenzahlfehler:  %.3e\n',N_relFehler);
fprintf('Max. Hermitizitätsfehler:      %.3e\n',hermFehler);
fprintf('Rel. Dichtefehler gg. exakt:   %.3e\n',n_relFehler);
fprintf('Rel. Stromfehler gg. exakt:    %.3e\n',j_relFehler);
% Teilchendichte
figure;
plot(X(:)*1e6, real(n_anfang(:)), 'y--', 'DisplayName','Anfangsdichte zu t=0'); hold on;
plot(X(:)*1e6, real(n_ende(:)), 'b.-', 'DisplayName','DG+FV');
plot(X(:)*1e6, real(n_ana(:)), 'r-', 'DisplayName','Analytisch');
xlabel('X [\mum]'); 
ylabel('n(X,\xi=0)');
legend('Location','best');
grid on;
title('Freier Dirac-Transport (beide Richtungen): Teilchendichte');
% Strom
figure;
plot(X(:)*1e6, real(j(:)), 'b.-', 'DisplayName','DG+FV'); hold on;
plot(X(:)*1e6, j_ana(:), 'r-', 'DisplayName','Analytisch');
xlabel('X [\mum]');
ylabel('j(X,t)');
legend('Location','best');
grid on;
title('Verlauf der Stromdichte');
% Gitterkonvergenz
figure;
loglog(h, n_Fehler, 'bo-', 'LineWidth', 1.5, 'DisplayName', 'Dichtefehler');
hold on;
loglog(h, j_Fehler, 'rs-', 'LineWidth', 1.5, 'DisplayName', 'Stromfehler');
xlabel('DG-Gitterweite \Delta X [m]');
ylabel('Relativer Fehler');
title('Gitterkonvergenz des DG-FV-Verfahrens');
legend('Location','best');
grid on;


