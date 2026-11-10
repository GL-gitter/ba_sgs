%% Schritt 0: freier 1D-Dirac-Transport einer 2x2-Dichtematrix
% DG in X, FV in xi, periodische Testgeometrie.
% Voraussetzung: JacobiGL.m, Vandermonde1D.m und Dmatrix1D.m im MATLAB-Pfad.
clear; clc;

%% 1) Gitter und DG-Matrizen
N = 1;                        % Zwei DG-Knoten pro Element
K = 12;                       % DG-Elemente in X
Nxi = 15;                     % FV-Zellen in xi (ungerade: xi=0 ist Zellmitte)
Lx = 4e-6;                    % Periodische Laenge in X [m]
dX = Lx / K;
d_xi = 2*Lx / Nxi;          % xi ist periodisch mit Periode 2*Lx

r = JacobiGL(0,0,N);
V = Vandermonde1D(N,r);
D = Dmatrix1D(N,r,V);
M = inv(V*V');
S = M*D;                     % Definition aus der Bachelorarbeit
M_k = (dX/2)*M;

% X-Knoten im jeweiligen DG-Element und xi-Zellmittelpunkte
X = zeros(2,K);
for k = 1:K
    X(:,k) = (k-1)*dX + (r+1)*dX/2;
end
xi = ((1:Nxi) - (Nxi+1)/2)*d_xi;

%% 2) Transportmatrizen fuer ky=0, Delta=0, B=0
vF = 1e6;                    % [m/s]
A = vF*[0 1; 1 0];
I = eye(2);
A_X = 0.5*(kron(I,A) + kron(A.',I));
A_xi = kron(I,A) - kron(A.',I);

[EV_X,EW_X] = eig(A_X);
lambda_X = diag(EW_X);
A_X_plus = EV_X*diag(max(lambda_X,0))*EV_X';
A_X_minus = EV_X*diag(min(lambda_X,0))*EV_X';

[EV_xi,EW_xi] = eig(A_xi);
lambda_xi = diag(EW_xi);
A_xi_plus = EV_xi*diag(max(lambda_xi,0))*EV_xi';
A_xi_minus = EV_xi*diag(min(lambda_xi,0))*EV_xi';

%% 3) Anfangszustand: Gauss-Wellenpaket mit beiden Laufrichtungen
% Psi(x,0) = g(x)*[1;0]. Dieser Zustand ist eine Ueberlagerung
% von rechts- und linkslaufenden Dirac-Moden: DG und FV sind aktiv.
% rho(X,xi,0) = Psi(X+xi/2)*Psi(X-xi/2)^dagger
x0 = 1.2e-6;                 % Position des Paketzentrums [m]
sigma = 0.4e-6;             % Paketbreite [m]
g = @(x) exp(-0.5*(mod(x-x0+Lx/2,Lx)-Lx/2).^2/sigma^2);
chi0 = [1;0];
rho_spin = chi0*chi0';

% U(:,n,k,j): 4 rho-Komponenten, DG-Knoten, X-Element, xi-Zelle
U = zeros(4,2,K,Nxi);
for k = 1:K
    for j = 1:Nxi
        for n = 1:2
            x1 = X(n,k) + xi(j)/2;
            x2 = X(n,k) - xi(j)/2;
            rho = g(x1)*g(x2)*rho_spin;
            U(:,n,k,j) = rho(:);    % [rho_AA; rho_BA; rho_AB; rho_BB]
        end
    end
end
U_anfang = U;

%% 4) Zeitentwicklung: explizites Runge-Kutta-Verfahren 4. Ordnung
T = 0.4e-12;                 % Endzeit [s]
dt_max = 0.05*min(dX/vF,d_xi/(2*vF));  % Konservativer Zeitschritt
nSteps = ceil(T/dt_max);
dt = T/nSteps;

rhs = @(W) DG_FV_rhs(W,S,M_k,d_xi, ...
    A_X,A_X_plus,A_X_minus,A_xi_plus,A_xi_minus);

for it = 1:nSteps
    k1 = rhs(U);
    k2 = rhs(U + 0.5*dt*k1);
    k3 = rhs(U + 0.5*dt*k2);
    k4 = rhs(U + dt*k3);
    U = U + (dt/6)*(k1+2*k2+2*k3+k4);
end

%% 5) Erste Kontrollen: Dichte, globale Bilanz, Hermitizitaet
j0 = (Nxi+1)/2;              % xi=0
n_anfang = squeeze(U_anfang(1,:,:,j0) + U_anfang(4,:,:,j0));
n_ende   = squeeze(U(1,:,:,j0) + U(4,:,:,j0));
mass_anfang = (dX/2)*sum(n_anfang(:));
mass_ende   = (dX/2)*sum(n_ende(:));
relMassfehler = abs(mass_ende-mass_anfang)/mass_anfang;

% Hermitizitaet: rho(X,xi) = rho(X,-xi)^dagger
hermFehler = 0;
for j = 1:Nxi
    for k = 1:K
        for n = 1:2
            rho_plus = reshape(U(:,n,k,j),2,2);
            rho_minus = reshape(U(:,n,k,Nxi+1-j),2,2);
            hermFehler = max(hermFehler,norm(rho_plus-rho_minus','fro'));
        end
    end
end

% Analytisch: anfangs spinpolarisiertes Paket trennt sich
% in eine rechts- und eine linkslaufende Dirac-Mode.
n_exakt = 0.5*(g(X-vF*T).^2 + g(X+vF*T).^2);
relDichtefehler = norm(n_ende(:)-n_exakt(:))/norm(n_exakt(:));

% Exakter Strom: rechtslaufende minus linkslaufende Mode
strom = squeeze(vF*(U(2,:,:,j0)+U(3,:,:,j0)));
strom_exakt = 0.5*vF*(g(X-vF*T).^2 - g(X+vF*T).^2);
relStromfehler = norm(strom(:)-strom_exakt(:))/norm(strom_exakt(:));

fprintf('Zeitschritte:                  %d\n',nSteps);
fprintf('Relativer Massenfehler:        %.3e\n',relMassfehler);
fprintf('Max. Hermitizitaetsfehler:     %.3e\n',hermFehler);
fprintf('Rel. Dichtefehler gg. exakt:   %.3e\n',relDichtefehler);
fprintf('Rel. Stromfehler gg. exakt:    %.3e\n',relStromfehler);

figure;
plot(X(:)*1e6, n_anfang(:), 'k--', 'DisplayName','Anfang'); hold on;
plot(X(:)*1e6, n_ende(:), 'b.-', 'DisplayName','DG+FV nach T');
plot(X(:)*1e6, n_exakt(:), 'r-', 'DisplayName','Analytisch nach T');
xlabel('X [um]'); ylabel('n(X,xi=0)');
legend('Location','best'); grid on;
title('Schritt 0: Freier Dirac-Transport (beide Richtungen)');

%% Lokale Funktion: rechte Seite fuer alle DG-Elemente und FV-Zellen
function dUdt = DG_FV_rhs(U,S,M_k,d_xi,A_X,A_X_plus,A_X_minus,A_xi_plus,A_xi_minus)
    K = size(U,3);
    Nxi = size(U,4);
    dUdt = zeros(size(U));

    for k = 1:K
        k_links = mod(k-2,K)+1;     % periodischer linker Nachbar
        k_rechts = mod(k,K)+1;     % periodischer rechter Nachbar

        for j = 1:Nxi
            j_links = mod(j-2,Nxi)+1;
            j_rechts = mod(j,Nxi)+1;

            U_mitte = U(:,:,k,j);  % 4x2-Matrix

            % DG-Randfluesse in X
            F_links = A_X_plus*U(:,2,k_links,j) ...
                    + A_X_minus*U_mitte(:,1);
            F_rechts = A_X_plus*U_mitte(:,2) ...
                     + A_X_minus*U(:,1,k_rechts,j);

            % Schwache DG-Formulierung wie in der Bachelorarbeit
            F = A_X*U_mitte;
            DG = (F*S - F_rechts*[0 1] + F_links*[1 0])/M_k;

            % FV-Flussdifferenz in xi
            G = (A_xi_plus*(U_mitte-U(:,:,k,j_links)) ...
               + A_xi_minus*(U(:,:,k,j_rechts)-U_mitte))/d_xi;

            % Zeitliche Ableitung: dU/dt = DG - G
            dUdt(:,:,k,j) = DG-G;
        end
    end
end
