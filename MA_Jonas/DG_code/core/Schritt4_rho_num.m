% Vorher "Schritt2.m" und "Schritt3.m" ausführen für die Variablen im
% Workspace
%% Definitionen
% Testgitter
N = 1; % DG-Polynomgrad
K = 8; % Ein DG-Element
N_xi = 17; % Drei FV-Zellen
d_X = L/K;
d_xi = 2*L/N_xi;


%% Physikalisches Gebiet und abgeschnittene Zellen
% Geometrie in nm fuer besser handhabbare Zahlen
L_nm = L*1e9;
X_kanten = linspace(0,L_nm,K+1);
xi_kanten = linspace(-L_nm,L_nm,N_xi+1);
% Rautenfoermiges physikalisches Gebiet
Gebiet = polyshape( ...
    [0, L_nm/2, L_nm, L_nm/2], ...
    [0, L_nm, 0, -L_nm]);
% Fuer jede Zelle den physikalisch zulaessigen Teil bestimmen
Zellen = cell(K,N_xi);
aktiv = false(K,N_xi);
for k = 1:K
    for j = 1:N_xi
        % Urspruengliche rechteckige DG/FV-Zelle
        Rechteck = polyshape( ...
            [X_kanten(k), X_kanten(k+1), ...
            X_kanten(k+1), X_kanten(k)], ...
            [xi_kanten(j), xi_kanten(j), ...
            xi_kanten(j+1), xi_kanten(j+1)]);
        % Schnittmenge mit dem physikalischen Gebiet
        Zellen{k,j} = intersect(Rechteck,Gebiet);
        % Hat die Zelle einen physikalisch zulaessigen Flaechenanteil?
        aktiv(k,j) = area(Zellen{k,j}) > 1e-10;
    end
end
% Physikalische Zellen darstellen
figure;
hold on;
for k = 1:K
    for j = 1:N_xi
        if aktiv(k,j)
            plot(Zellen{k,j});
        end
    end
end
axis equal;
xlabel('X [nm]');
ylabel('\xi [nm]');
title('Physikalisches DG/FV-Gitter');

%% Aktive Zellen nummerieren
Zellennummer = zeros(K,N_xi);
Zellennummer(aktiv) = 1:nnz(aktiv);
% 16 Komponenten pro DG-Knoten,
% N+1 DG-Knoten pro aktiver Zelle
N_u = 16*(N+1)*nnz(aktiv);

K_Spektral_1 = sparse(N_u,N_u);
K_Spektral_2 = sparse(N_u,N_u);

r = JacobiGL(0,0,N);
V = Vandermonde1D(N,r);
D = Dmatrix1D(N,r,V);
M = inv(V*V');
S = M*D;
M_k = (d_X/2)*M;
Z = zeros(8,16);
I4 = eye(4);
I2 = eye(2);
K_mitte_alle = cell(2,length(phi_werte));

s_VZ = diag([1,-1]);
A = vF*kron(s_VZ,sigma_x);
% Aufspaltung A = A_plus + A_minus
[W,Lambda] = eig(A);
lambda = real(diag(Lambda));
A_plus = W*diag(max(lambda,0))*W';
A_minus = W*diag(min(lambda,0))*W';

for i=1:length(phi_werte)
    phi = phi_werte(i);
    for n=1:2
        E_n = E_ABS_plus_num(i,n);
        ML = M_L(E_n,phi);
        MR = M_R(E_n,phi);


        % Andreev-Randbedingung fuer das erste Ortsargument x1
        R1_L = kron(I4,[-ML, I2]); % Linker Kontakt x1=0
        R1_R = kron(I4,[-MR, I2]); % Rechter Kontakt x1=L
        % Andreev-Randbedingung fuer das zweite Ortsargument x2
        Q_L = [-ML'; I2];
        Q_R = [-MR'; I2];
        R2_L = kron(Q_L.',I4);  % Linker Kontakt x2=0
        R2_R = kron(Q_R.',I4);  % Rechter Kontakt x2=L
        % Andreev-Randbedingungen auf einem 2x2-Gitter in (x1,x2)
        K_RB = [ ...
            R1_L, Z, Z, Z; ...  % x1=0, x2=0
            Z, R1_R, Z, Z; ...  % x1=L, x2=0
            Z, Z, R1_L, Z; ...  % x1=0, x2=L
            Z, Z, Z, R1_R; ...  % x1=L, x2=L

            R2_L, Z, Z, Z; ...  % x2=0, x1=0
            Z, R2_L, Z, Z; ...  % x2=0, x1=L
            Z, Z, R2_R, Z; ...  % x2=L, x1=0
            Z, Z, Z, R2_R  ...  % x2=L, x1=L
            ];
        %size(K_RB)  % 64 x 64


        % 1. Term, Vektorisierung von A*rho
        S_weak = S.';
        C = kron(I4,A);
        DG_1 = (1i*hbar/2)*[S_weak(1,1)*C, S_weak(1,2)*C];
        DG_2 = (1i*hbar/2)*[S_weak(2,1)*C, S_weak(2,2)*C];
        % 3. Term, FV-Flüsse
        Cp = kron(I4,A_plus);
        Cm = kron(I4,A_minus);
        FV = [-Cp, Cp-Cm, Cm];
        % Hamilton-Potential im Graphen (ky = 0)
        B = kron(s_VZ, (U0-mu)*I2);
        % Potentialmatrizen für i=1
        P11 = M_k(1,1)*(B-E_n*I4);
        P12 = M_k(1,2)*(B-E_n*I4);
        % Potentialmatrizen für i=2
        P21 = M_k(2,1)*(B-E_n*I4);
        P22 = M_k(2,2)*(B-E_n*I4);

        Pot_1 = [kron(I4,P11), kron(I4,P12)];
        Pot_2 = [kron(I4,P21), kron(I4,P22)];
        K_Pot_1 = [zeros(16,32), Pot_1, zeros(16,32)];
        K_Pot_2 = [zeros(16,32), Pot_2, zeros(16,32)];

        % DG-Randflussterm
        % Testknoten i=1:
        % +ihbar/2 * A_minus * rho_(1,2,1)
        K_Rand_1 = (1i*hbar/2)*[ ...
            zeros(16,32), Cm, zeros(16,16), zeros(16,32)];
        % Testknoten i=2:
        % -ihbar/2 * A_plus * rho_(1,2,2)
        K_Rand_2 = -(1i*hbar/2)*[ ...
            zeros(16,32), zeros(16,16), Cp, zeros(16,32)];

        % Erster Testknoten i=1
        K_DG_1 = [zeros(16,32), DG_1, zeros(16,32)];
        K_FV_1 = -(1i*hbar/d_xi)*[ ...
            -M_k(1,1)*Cp, ...
            -M_k(1,2)*Cp, ...
            M_k(1,1)*(Cp-Cm), ...
            M_k(1,2)*(Cp-Cm), ...
            M_k(1,1)*Cm, ...
            M_k(1,2)*Cm  ...
            ];

        % Zweiter Testknoten i=2
        K_DG_2 = [zeros(16,32), DG_2, zeros(16,32)];
        K_FV_2 = -(1i*hbar/d_xi)*[ ...
            -M_k(2,1)*Cp, ...
            -M_k(2,2)*Cp, ...
            M_k(2,1)*(Cp-Cm), ...
            M_k(2,2)*(Cp-Cm), ...
            M_k(2,1)*Cm, ...
            M_k(2,2)*Cm  ...
            ];

        K_1 = K_DG_1 + K_Rand_1 + K_FV_1 + K_Pot_1;
        K_2 = K_DG_2 + K_Rand_2 + K_FV_2 + K_Pot_2;
        K_mitte = [K_1; K_2];
        K_mitte_alle{n,i} = K_mitte;
        K_gesamt = 0;
    end
end