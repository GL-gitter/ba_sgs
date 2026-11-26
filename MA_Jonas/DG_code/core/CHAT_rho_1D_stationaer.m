clc, clearvars
%% Matrizen und Parameter
N = 1;                      % DG-Polynomgrad (hier fest N=1)
K = 12;                     % Anzahl DG-Elemente
L_X = 4e-6;                 % Teststrecke [m]
d_X = L_X/K;
vF = 1e6;                   % Fermi-Geschwindigkeit [m/s]

r = JacobiGL(0,0,N);
V = Vandermonde1D(N,r);
D = Dmatrix1D(N,r,V);
M = inv(V*V');
S = M*D;
M_k = (d_X/2)*M;

X = zeros(N+1,K);
for k = 1:K
    X(:,k) = 0.5*d_X*r + (k-0.5)*d_X;
end

N_kx = 81;
kx_werte = linspace(-2e7,2e7,N_kx);  % [1/m]

%% Dirac-Transportmatrizen (wie in Schritt0.m)
A = vF*[0 1; 1 0];
I = eye(2);
A_X = 0.5*(kron(I,A) + kron(A.',I));
A_xi = kron(I,A) - kron(A.',I);

[EV_X,EW_X] = eig(A_X);
lambda_X = diag(EW_X);
lambda_X(abs(lambda_X) < 1e-12*vF) = 0;
A_X_plus = EV_X*diag(max(lambda_X,0))*EV_X';
A_X_minus = EV_X*diag(min(lambda_X,0))*EV_X';
P0 = EV_X(:,lambda_X==0)*EV_X(:,lambda_X==0)';

%% Einlaufende Randverteilungen (frei gewaehlter Test)
chi_plus = [1;1]/sqrt(2);
chi_minus = [1;-1]/sqrt(2);
P_plus = chi_plus*chi_plus';
P_minus = chi_minus*chi_minus';

k0 = 1e7;                   % [1/m]
sigma_k = 2e6;              % [1/m]
f_L = @(kx) exp(-0.5*((kx-k0)/sigma_k).^2);
f_R = @(kx) 0.5*exp(-0.5*((kx+k0)/sigma_k).^2);

%% DG-Systemmatrix (schwache Form, offene Upwind-Randfluesse)
e_L = [1;0];
e_R = [0;1];
nd = 4*(N+1);              % Freiheitsgrade pro Element = 8
nd_ges = nd*K;

K_mitte = -kron(S.',A_X) ...
    + kron(e_R*e_R.',A_X_plus) ...
    - kron(e_L*e_L.',A_X_minus);
K_links = -kron(e_L*e_R.',A_X_plus);
K_rechts = kron(e_R*e_L.',A_X_minus);

K_DG = sparse(nd_ges,nd_ges);
for k = 1:K
    idx = (k-1)*nd+(1:nd);
    K_DG(idx,idx) = K_mitte;
    if k > 1
        K_DG(idx,idx-nd) = K_links;
    end
    if k < K
        K_DG(idx,idx+nd) = K_rechts;
    end
end
K_xi = kron(speye(K),sparse(kron(M_k.',A_xi)));

%% Stationaere Loesung fuer jeden kx-Wert
W = complex(zeros(4,N+1,K,N_kx));
f = complex(zeros(2,2,N+1,K,N_kx));
Fehler = zeros(1,N_kx);
Res = zeros(1,N_kx);
herm_Fehler = zeros(1,N_kx);

for j = 1:N_kx
    kx = kx_werte(j);
    f_ein_L = f_L(kx)*P_plus;
    f_ein_R = f_R(kx)*P_minus;
    w_L = f_ein_L(:);
    w_R = f_ein_R(:);

    b = zeros(nd_ges,1);
    b(1:nd) = kron(e_L,A_X_plus*w_L);
    b(end-nd+1:end) = b(end-nd+1:end) ...
        - kron(e_R,A_X_minus*w_R);

    K_phys = K_DG + 1i*kx*K_xi;
    K_sys = K_phys;
    if kx == 0
        % Nullgeschwindigkeits-Kohaerenzen sind bei kx=0 unbestimmt.
        % Zusatzannahme fuer den inkohärenten Test: P0*w=0 an allen Knoten.
        K_sys = K_sys + vF*kron(speye((N+1)*K),sparse(P0));
    end

    w = (K_sys/vF)\(b/vF);
    W(:,:,:,j) = reshape(w,4,N+1,K);
    f(:,:,:,:,j) = reshape(w,2,2,N+1,K);

    % Analytische Vergleichsloesung: konstanter, ballistischer Transport
    f_ana = f_ein_L + f_ein_R;
    w_ana = repmat(f_ana(:), (N+1)*K,1);
    Fehler(j) = norm(w-w_ana)/max(norm(w_ana),eps);
    Res(j) = norm(K_phys*w-b)/max(norm(b),eps);
    for k = 1:K
        for a = 1:N+1
            F = f(:,:,a,k,j);
            herm_Fehler(j) = max(herm_Fehler(j),norm(F-F','fro'));
        end
    end
end

%% Rueckgewinnung von rho(X,xi=0), Teilchen- und Stromdichte
% Fourier-Konvention: f = (1/2pi)*integral rho*exp(-ikx*xi) dxi.
% Damit rho(X,0) = integral f(X,kx) dkx.
rho = trapz(kx_werte,f,5);
n = zeros(N+1,K);
j_strom = zeros(N+1,K);
for k = 1:K
    for a = 1:N+1
        R = rho(:,:,a,k);
        n(a,k) = real(trace(R));
        j_strom(a,k) = real(trace(A*R));
    end
end

n_ana = trapz(kx_werte,f_L(kx_werte)+f_R(kx_werte));
j_ana = vF*trapz(kx_werte,f_L(kx_werte)-f_R(kx_werte));

fprintf('Max. relativer Fehler:       %.3e\n',max(Fehler));
fprintf('Max. relatives Residuum:     %.3e\n',max(Res));
fprintf('Max. Hermitizitaetsfehler:    %.3e\n',max(herm_Fehler));
fprintf('Relativer Dichtefehler:      %.3e\n',max(abs(n(:)-n_ana))/abs(n_ana));
fprintf('Relativer Stromfehler:       %.3e\n',max(abs(j_strom(:)-j_ana))/abs(j_ana));
