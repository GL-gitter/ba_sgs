function [G_glob, info] = get_Drift_rho(mat, p, U_IV, EM)
%GET_DRIFT_RHO Potential/drift operator assembled in rho basis.
%
% This keeps the potential action in the relative-coordinate basis instead
% of transforming it into the phase/eigen basis with phi'*C*phi.

if ~p.permute_shell
    error('get_Drift_rho requires mat.dg.params.permute_shell = true.');
end

if (p.doIV == true && p.doPoisson == false)
    p.U = U_IV;
    p.U1 = zeros(size(p.xg));
    idxRamp = p.xg < p.L_chi/2+p.L_d/2 & p.xg > p.L_chi/2-p.L_d/2;
    p.U1(idxRamp) = linspace(0, U_IV, sum(idxRamp));
    p.U1(p.xg >= p.L_chi/2+p.L_d/2) = U_IV;
end
if (p.doPoisson == true && p.doIV == true)
    p.U = U_IV;
end

p.M = p.delta_xi/2*p.M;
p.invM = inv(p.M);

x_anteil = diff(p.x);
x_prozent = x_anteil/2;
chi_anteil = x_prozent*p.delta_chi;
xi_anteil = x_prozent*p.delta_xi;
chi_GL_koord = zeros(1, p.Nqx);
xi_GL_koord = zeros(1, p.Nqx);
for i = 1:p.Nqx-1
    chi_GL_koord(i+1) = chi_GL_koord(i) + chi_anteil(i);
    xi_GL_koord(i+1) = xi_GL_koord(i) + xi_anteil(i);
end

if (p.doFV == true)
    xiKantenkoord_GL = linspace(-p.L_xi/2, p.L_xi/2, p.N_xi+1);
elseif (p.doDG == true)
    xiKantenkoord_GL = zeros(1, p.Ny);
    for i = 0:p.N_xi-1
        xiKantenkoord_GL(1, i*p.Nqx+1:(i+1)*p.Nqx) = p.xi_Rechengebiet(i+1) + xi_GL_koord;
    end
elseif (p.doFEM == true)
    xiKantenkoord_GL = zeros(1, p.Ny);
    for i = 0:p.N_xi
        xiKantenkoord_GL(1, i*p.Nqx+1:(i+1)*p.Nqx) = p.xi_Rechengebiet(i+1) + xi_GL_koord;
    end
else
    error('Unknown xi discretization.');
end

y_drift = zeros(p.Nqx*p.N_chi, p.Ny);
for i = 1:p.Nqx*p.N_chi
    y_drift(i, :) = xiKantenkoord_GL;
end

Kantenkoord_GL = zeros(1, p.N_chi*p.Nqx);
for i = 0:p.N_chi-1
    Kantenkoord_GL(1, i*p.Nqx+1:(i+1)*p.Nqx) = p.disk_Rechengebiet1(i+1) + chi_GL_koord;
end

x_drift = zeros(p.N_chi*p.Nqx, p.Ny);
for i = 1:p.Ny
    x_drift(:, i) = Kantenkoord_GL.';
end

Bpotential = V_transient(x_drift + 0.5*y_drift, mat.x*1e-9, EM, p.U, p.L_chi) ...
           - V_transient(x_drift - 0.5*y_drift, mat.x*1e-9, EM, p.U, p.L_chi);
Bcap = get_CAP(y_drift, p.L_xi, p.delta_xi);
B = Bpotential - 1i*Bcap;

Nrel = size(p.phi, 1);
W = p.V_tilde/p.V;

if (p.doFV == true)
    maxPairs = max(1, 3*Nrel - 2);
elseif (p.doDG == true)
    maxPairs = max(1, p.N_xi*p.N_K_xi*p.N_K_xi);
else
    maxPairs = Nrel*Nrel;
end

maxEntries = p.N_chi*maxPairs*p.N_K_chi*p.N_K_chi;
rowIdx = zeros(maxEntries, 1);
colIdx = zeros(maxEntries, 1);
values = complex(zeros(maxEntries, 1));
entry = 0;

[rowLocal, colLocal] = ndgrid(1:p.N_K_chi, 1:p.N_K_chi);
tol = max(1e-15, 1e-16*norm(B(:)));

for d = 1:p.N_chi
    Gd = zeros(p.Nqx, Nrel, Nrel);
    for q = 1:p.Nqx
        kk = (d-1)*p.Nqx + q;
        Gd(q, :, :) = localPotentialRho(B(kk, :), p, Nrel);
    end

    pairMask = squeeze(any(abs(Gd) > tol, 1));
    [relRows, relCols] = find(pairMask);

    for pairIdx = 1:length(relRows)
        j = relRows(pairIdx);
        m = relCols(pairIdx);
        c = squeeze(Gd(:, j, m)).*p.w;
        if norm(c) <= tol
            continue
        end

        W_tilde = W.*c;
        Drift = p.Q_drift*p.invM.*(W_tilde.'*W);

        rows = localToGlobal(d, rowLocal, j, p.N_K_chi, Nrel);
        cols = localToGlobal(d, colLocal, m, p.N_K_chi, Nrel);
        nAdd = numel(Drift);
        if entry + nAdd > numel(values)
            growBy = max(maxEntries, nAdd);
            rowIdx = [rowIdx; zeros(growBy, 1)]; %#ok<AGROW>
            colIdx = [colIdx; zeros(growBy, 1)]; %#ok<AGROW>
            values = [values; complex(zeros(growBy, 1))]; %#ok<AGROW>
        end

        range = entry + (1:nAdd);
        rowIdx(range) = rows(:);
        colIdx(range) = cols(:);
        values(range) = Drift(:);
        entry = entry + nAdd;
    end
end

rowIdx = rowIdx(1:entry);
colIdx = colIdx(1:entry);
values = values(1:entry);
G_glob = sparse(rowIdx, colIdx, values, p.N_chi*p.N_K_chi*Nrel, p.N_chi*p.N_K_chi*Nrel);

info = struct;
info.Nrel = Nrel;
info.nnz = nnz(G_glob);
info.potentialChiNorm = zeros(p.N_chi, 1);
info.capChiNorm = zeros(p.N_chi, 1);
for d = 1:p.N_chi
    rows = (d-1)*p.Nqx+1:d*p.Nqx;
    info.potentialChiNorm(d) = norm(Bpotential(rows, :), 'fro');
    info.capChiNorm(d) = norm(Bcap(rows, :), 'fro');
end
info.potentialChiNormRange = [min(info.potentialChiNorm), max(info.potentialChiNorm)];
info.capChiNormRange = [min(info.capChiNorm), max(info.capChiNorm)];
end

function C = localPotentialRho(Brow, p, Nrel)
if (p.doFV == true)
    C = diag(Brow(1:end-1) + Brow(2:end)) ...
      + diag(Brow(2:end-1), -1) ...
      + diag(Brow(2:end-1), +1);
    C = 1i*(p.delta_xi/4)*C;
elseif (p.doDG == true)
    W_xi = p.V_tilde_xi/p.V_xi;
    W_xiCell = repmat({W_xi}, 1, p.N_xi);
    W_xi = blkdiag(W_xiCell{:});
    w_xi = repmat(p.w, p.N_xi, 1);
    c_xi = Brow(:).*w_xi;
    W_xi_tilde = W_xi.*c_xi;
    C = (1/1i)*p.M_xi_glob*W_xi_tilde.'*W_xi*p.delta_xi;
else
    error('get_Drift_rho currently supports FV and DG in xi.');
end

if ~isequal(size(C), [Nrel, Nrel])
    error('rho potential block has size [%d, %d], expected [%d, %d].', size(C, 1), size(C, 2), Nrel, Nrel);
end
end

function idx = localToGlobal(elem, localNode, relNode, N_K_chi, Nrel)
idx = (elem-1)*N_K_chi*Nrel + (localNode-1)*Nrel + relNode;
end
