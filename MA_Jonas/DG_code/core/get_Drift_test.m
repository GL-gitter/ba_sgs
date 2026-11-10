function G_glob = get_Drift_test(p,U_IV)
if (p.doIV == true)
    p.U = U_IV;
    p.U1 = zeros(size(p.xg));
    p.U1(p.xg < p.L_chi/2+p.L_d/2 & p.xg > p.L_chi/2-p.L_d/2) = linspace(0, p.U, sum(p.xg < p.L_chi/2+p.L_d/2 & p.xg > p.L_chi/2-p.L_d/2));
    p.U1(p.xg >= p.L_chi/2+p.L_d/2) = p.U;
end
p.M = p.delta_xi/2*p.M;
p.invM = inv(p.M);

%% Gauß-Lobatto coordinates
x_anteil = diff(p.x);                 % Abstände zwischen den GL-Punkten
x_prozent = x_anteil/2;             % Prozentualer Anteil der Quadraturpunkte
chi_anteil = x_prozent * p.delta_chi;
chi_GL_koord = zeros(1, p.Nqx);
xi_GL_koord = zeros(1, p.Nqx);
xi_anteil = x_prozent * p.delta_xi;
% Quadratur-Koordinaten für ein Xi-Element bestimmen
for i = 1:p.Nqx-1
    xi_GL_koord(i+1) = xi_GL_koord(i)+xi_anteil(i);
end
% Quadratur-Koordinaten für alle Xi-Elemente bestimmen
if (p.doFV == true)
    xiKantenkoord_GL = linspace(-p.L_xi/2, p.L_xi/2, p.N_xi+1);
elseif (p.doDG == true)
    xiKantenkoord_GL = zeros(1, p.Ny);
    for i = 0 : p.N_xi-1
        xiKantenkoord_GL(1, i*p.Nqx+1:(i+1)*p.Nqx) = p.xi_Rechengebiet(i+1)+xi_GL_koord;
    end
elseif (p.doFEM == true)
    xiKantenkoord_GL = zeros(1, p.Ny);
    for i = 0 : p.N_xi
        xiKantenkoord_GL(1, i*p.Nqx+1:(i+1)*p.Nqx) = p.xi_Rechengebiet(i+1)+xi_GL_koord;
    end
end
% xi-Koordinaten für den Driftoperator
y_drift = zeros(p.N_chi*p.Nqx, p.Ny);
for i = 1 : p.Nqx*p.N_chi
    y_drift(i, :) = xiKantenkoord_GL;
end

% Quadratur-Koordinaten für ein Chi-Element bestimmen
for i = 1:p.Nqx-1
    chi_GL_koord(i+1) = chi_GL_koord(i)+chi_anteil(i);
end

% Quadratur-Koordinaten für alle Chi-Elemente bestimmen
Kantenkoord_GL = zeros(1, p.N_chi*p.Nqx);
for i = 0 : p.N_chi-1
    Kantenkoord_GL(1, i*p.Nqx+1 : (i+1)*p.Nqx) = p.disk_Rechengebiet1(i+1) + chi_GL_koord;
end

% chi-Koordinaten für den Driftoperator
x_drift = zeros(p.N_chi*p.Nqx, p.Ny);
for i = 1 : p.Ny
    x_drift(:, i) = Kantenkoord_GL.';
end

% Driftwerte bestimmen und in einer Matrix (p.Nqx*p.N_chi, p.N_xi, p.N_xi) aufnehmen
B = Potential(x_drift+0.5*y_drift, p) - V_transient(x_drift+0.5*y_drift, p.xg, p.U1, p.U, p.L_chi) ...
    - Potential(x_drift-0.5*y_drift, p) + V_transient(x_drift-0.5*y_drift, p.xg, p.U1, p.U, p.L_chi) ...
    - 1i*get_CAP(y_drift, p.L_xi, p.delta_xi);

if (p.doDG == true)
    W_xi = p.V_tilde_xi/p.V_xi;
    W_xiCell = repmat({W_xi},1,p.N_xi);
    W_xi = blkdiag(W_xiCell{:});
    w_xi = repmat(p.w,p.N_xi,1);
    G = zeros(p.Nqx*p.N_chi, p.N_xi*p.N_K_xi, p.N_xi*p.N_K_xi);
    for kk = 1 : p.Nqx*p.N_chi
        c_xi = B(kk,:)'.*w_xi;
        W_xi_tilde = W_xi.*c_xi;
        C = p.M_xi_glob*W_xi_tilde.'*W_xi*p.delta_xi;
        G(kk,:,:) = 1/1i*p.phi'*C*p.phi;
    end
elseif (p.doFV == true)
    G = zeros(p.Nqx*p.N_chi, p.N_xi, p.N_xi);
    for i = 1 : p.Nqx*p.N_chi
        C = diag(B(i,1:end-1) + B(i,2:end))+diag(B(i,2:end-1),-1)+ diag(B(i,2:end-1),+1);
        C = p.delta_xi/4*C;
        G(i,:,:) = 1i*p.phi'*C*p.phi;
    end
elseif (p.doFEM == true)
    W_xi = p.V_tilde_xi/p.V_xi;
    W_xiCell = repmat({W_xi},1,p.N_xi+1);
    W_xi = blkdiag(W_xiCell{:});
    w_xi = repmat(p.w,p.N_xi+1,1);
    G = zeros(p.Nqx*p.N_chi, p.N_xi*(p.N_K_xi-1)+1, p.N_xi*(p.N_K_xi-1)+1);
    % Driftwerte in einer Tri-diagonalen Matrix aufnehmen und diagonalisieren
    for ii = 1 : p.Nqx*p.N_chi
        C_xi = zeros(p.N_xi*(p.N_K_xi-1)+1,p.N_xi*(p.N_K_xi-1)+1);
        c_xi = B(ii,:)'.*w_xi;
        W_xi_tilde = W_xi.*c_xi;
        C = W_xi_tilde.'*W_xi;
        C_xi(1:p.N_K_xi,1:p.N_K_xi)=C(1:p.N_K_xi,1:p.N_K_xi);
        for aa = 2:p.N_xi-1
            C_xi(aa:aa+p.N_K_xi-1, aa:aa+p.N_K_xi-1)=[C_xi(aa,aa)+C((aa-1)*p.N_K_xi+1,(aa-1)*p.N_K_xi+1), C((aa-1)*p.N_K_xi+1,(aa-1)*p.N_K_xi+p.N_K_xi);
                C((aa-1)*p.N_K_xi+p.N_K_xi,(aa-1)*p.N_K_xi+1),      C((aa-1)*p.N_K_xi+p.N_K_xi,(aa-1)*p.N_K_xi+p.N_K_xi)];
        end
        C_FEM = p.delta_xi/4*p.M_xi_glob*C_xi;
        G(ii,:,:) = 1/1i*p.phi'*C_FEM*p.phi;
    end
end

W = p.V_tilde/p.V;
W_all = repmat(W,p.N_chi,p.N_chi);
c = zeros(p.N_chi*p.Nqx, 1);
index = (1:p.N_K_chi*p.N_K_chi)';
index_zl = zeros(p.dim*p.dim/p.N_chi,1);
index_sp = index_zl;
werte_Drift = index_zl;
tol = 1e-16*norm(G(:));

if (p.doDG == true)
    Ny = p.N_xi*p.N_K_xi;
elseif (p.doFEM == true)
    Ny = p.N_xi*(p.N_K_xi-1)+1;
elseif (p.doFV == true)
    Ny = p.N_xi;
end
% w_all = repmat(p.w,p.N_chi,1);
% invM_allCell = repmat({p.invM},1,p.N_chi);
% invM_all = blkdiag(invM_allCell{:});
% for j = 1 : Ny
%     zl = ((j-1)*(p.N_K_chi*p.N_chi)+1 : j*(p.N_K_chi*p.N_chi))';
%     for m = 1 : Ny
%         sp = (m-1)*(p.N_K_chi*p.N_chi)+1 : m*(p.N_K_chi*p.N_chi);
%         % for d = 1 : p.N_chi
%             c(:, 1) = G(:, j, m).*w_all;      % Gauß-Lobatto Quadratur
%             if (norm(c) > tol)          % nummerische Nullen ausgrenzen
%                 xBlock = 1 : p.N_chi*p.N_K_chi;
%                 zl_inner = zl(xBlock)*ones(1, p.N_K_chi*p.N_chi);
%                 sp_inner = ones(p.N_K_chi*p.N_chi, 1)*(sp(xBlock));
%                 W_tilde = W_all.*c;
%                 Drift = p.Q_drift*invM_all.*(W_tilde.'*W_all);
%                 werte_Drift(index) = Drift(:);
%                 index_zl(index) = zl_inner;
%                 index_sp(index) = sp_inner;
%                 index = index + p.N_K_chi^2*p.N_chi^2;
%             end
%         % end
%     end
% end
% G_glob = [index_zl(:), index_sp(:), werte_Drift(:)];
% G_glob = myspconvert(G_glob, p.dim, p.dim, 1e-15);       % Aufstellen der glob. Driftmatrix und Ignorieren nummerischer Nullen (1e-15)

G_glob = zeros(p.dim,p.dim);
G_j = zeros(Ny,p.N_K_chi*p.N_chi,p.dim);
parfor j = 1:Ny
    G_j(j,:,:) = driftmatrix_loc(p,j,werte_Drift,index_zl,index_sp,index,W,tol,Ny,G);
end
for j = 1:Ny
    G_glob((j-1)*p.N_K_chi*p.N_chi+1:j*p.N_K_chi*p.N_chi,:) = G_j(j,:,:);
end
G_glob = sparse(G_glob);
end
