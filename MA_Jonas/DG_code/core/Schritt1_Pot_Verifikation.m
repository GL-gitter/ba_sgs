% Erst "Schritt1.m" ausführen für Variable U (mit Potential) im Workspace

%% Matrizen und Parameter
N = 1; % Polynomgrad des DG-Ansatzes, N+1 Knoten
K_werte = [12, 24, 48]; %, 96]; % Unterteilung der X-Richtung in 12 DG-Elemente
N_xi_werte = [15, 31, 61]; %, 121]; % Unterteilung der xi-Richtung in 15 FV-
                                % Zellen, ungerade für Mittelwert bei 0                               
U0 = 0; % Barrierenhöhe: 2 meV [J]

U_mit_Pot = U;
U_ohne_Pot = @(x) U0*exp(-0.5*((x-x_B)/w_B).^2);

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
    t_End = 1.8e-12; % Festlegung der Endzeit (0.4ps für 0.4µm) für RK4 später
    
    
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
    S_pot_ohne = zeros(2,K,N_xi);


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
    
    %% Zeitentwicklung: explizites Runge-Kutta-Verfahren 4. Ordnung
    dt_max = 0.05*min(d_X/vF,d_xi/(2*vF)); % CFL Bedingung maximal möglicher 
                % Zeitschritt, c=0.05 (Testwert). Minimum beider Werte, damit 
                % keine Zelle ggf. übersprungen wird.
    n_steps = ceil(t_End/dt_max); % Anzahl der Zeitschritte insgesamt
    dt = t_End/n_steps; % Schrittweite (h_k im Skript)
    % Definition des zu lösenden Problems und Lösungsschleife
    rechts_ohne = @(W) DG_FV_rechts(W, S, M_k, d_xi, A_X, A_X_plus, A_X_minus, ...
        A_xi_plus, A_xi_minus, S_pot_ohne);
    for it = 1:n_steps
        k1 = rechts_ohne(U);
        k2 = rechts_ohne(U + 0.5*dt*k1);
        k3 = rechts_ohne(U + 0.5*dt*k2);
        k4 = rechts_ohne(U + dt*k3);
        U = U + (dt/6)*(k1+2*k2+2*k3+k4);
    end
end
U_ohne_Potential = U;

%% Ausgabe
differenz = U_mit_Pot - U_ohne_Potential;

rel_Unterschied = norm(differenz(:))/norm(U_ohne_Potential(:));

fprintf('Relativer Unterschied durch Gaußbarriere: %.3e\n', ...
    rel_Unterschied);