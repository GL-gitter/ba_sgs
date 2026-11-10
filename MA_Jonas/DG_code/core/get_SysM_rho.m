function [A, rhs, info] = get_SysM_rho(mat, p, EM)
%GET_SYSM_RHO System matrix for the rho-basis DG transport variant.

[A, rhs, fluxInfo] = get_Diff_rho(mat, p);
[G_glob, driftInfo] = get_Drift_rho(mat, p, 0, EM);

diffScale = 1;
driftScale = 1;
if isfield(mat, 'dg') && isfield(mat.dg, 'params')
    if isfield(mat.dg.params, 'rho_diff_scale')
        diffScale = mat.dg.params.rho_diff_scale;
    end
    if isfield(mat.dg.params, 'rho_drift_scale')
        driftScale = mat.dg.params.rho_drift_scale;
    end
end

A_diff = diffScale*A;
G_glob = driftScale*G_glob;
A = A_diff + G_glob;

info = struct;
info.flux = fluxInfo;
info.drift = driftInfo;
info.scale.diff = diffScale;
info.scale.drift = driftScale;
% info.nnzDiffAndFlux = nnz(A_diff);
% info.nnzDrift = nnz(G_glob);
% info.normestDiffAndFlux = normest(A_diff);
% info.normestDrift = normest(G_glob);
% info.normestDriftToDiff = info.normestDrift/max(info.normestDiffAndFlux, eps);
end
