function [A, rhs, fluxInfo] = get_Diff_rho(mat, p)
%GET_DIFF_RHO DG transport operator in chi with rho as global unknown.
%
% The old solver stores the global unknown in the diagonal phase/eigen basis.
% This variant stores blocks directly in the relative-coordinate basis:
%
%   rho_global = [rho(chi_1); rho(chi_2); ...]
%
% fluxType = 'central' uses the sparse local rho operator -1i*d/dxi.
% The legacy dense characteristic matrix-upwind path is kept separate.
% The dense projected operator Phi*D*PhiLeft is available separatly as
% 'central-projected' for diagnostics. The domain boundaries still use the
% characteristic inflow contribution.
%%



if ~p.permute_shell
    error('get_Diff_rho requires mat.dg.params.permute_shell = true.');
end

fluxType = 'central';
if isfield(mat, 'dg') && isfield(mat.dg, 'params')
    if isfield(mat.dg.params, 'rho_flux')
        fluxType = mat.dg.params.rho_flux;
    end
end

fnl_c = p.verteilung_l(:);
fnr_c = p.verteilung_r(:);

Dr = p.M*Dmatrix1D(p.N_K_chi-1, p.r, p.V);
detJx = ones(p.N_chi, 1)*p.delta_chi/2;
p.M = p.delta_xi/2*p.M;
p.invM = inv(p.M);
Dr(abs(Dr) < 1e-15) = 0;

P1 = zeros(p.N_K_chi);
P1(1, 1) = 1;
P2 = zeros(p.N_K_chi);
P2(1, p.N_K_chi) = 1;
P3 = zeros(p.N_K_chi);
P3(p.N_K_chi, p.N_K_chi) = 1;
P4 = zeros(p.N_K_chi);
P4(p.N_K_chi, 1) = 1;

K1 = p.invM*(Dr + 0.5*P1 - 0.5*P3);
K2 = p.invM*(0.5*P1 + 0.5*P3);
K2L = p.invM*(0.5*P1);
K2R = p.invM*(0.5*P3);
K3 = -0.5*p.invM*P2;
K4 = K3;
K5 = 0.5*p.invM*P4;
K6 = -K5;

H1 = kron(spdiags(1./detJx, 1, p.N_chi, p.N_chi), K5) ...
   + kron(spdiags(1./detJx(2:end), -1, p.N_chi, p.N_chi), K3) ...
   + kron(spdiags(1./detJx, 0, p.N_chi, p.N_chi), K1);
H2 = kron(spdiags(1./detJx, 1, p.N_chi, p.N_chi), K6) ...
   + kron(spdiags(1./detJx(2:end), -1, p.N_chi, p.N_chi), K4) ...
   + kron(spdiags(1./detJx, 0, p.N_chi, p.N_chi), K2);
H2Boundary = boundaryUpwindChiMatrix(p, detJx, K2L, K2R);

Phi = p.phi;
PhiLeft = leftInversePhi(Phi);
lambda = diag(p.D);
lambdaPlus = (real(lambda) > 0).*lambda;
lambdaMinus = (real(lambda) < 0).*lambda;

LambdaPlus = spdiags(lambdaPlus, 0, length(lambda), length(lambda));
LambdaMinus = spdiags(lambdaMinus, 0, length(lambda), length(lambda));

Aplus = Phi*LambdaPlus*PhiLeft;
Aminus = Phi*LambdaMinus*PhiLeft;
Achi = Aplus + Aminus;
Aabs = Aplus -Aminus;

rhoL = Phi*fnl_c;
rhoR = Phi*fnr_c;

Q1 = zeros(p.N_K_chi*p.N_chi, 1);
Q1(1:p.N_K_chi) = p.invM(:, 1)/detJx(1);
Q2 = zeros(p.N_K_chi*p.N_chi, 1);
Q2(p.N_K_chi*(p.N_chi-1)+1:p.N_K_chi*p.N_chi) = p.invM(:, p.N_K_chi)/detJx(end);

switch lower(fluxType)
    case {'matrix-upwind', 'upwind', 'characteristic-upwind'}
        A = p.Q_diff*(kron(H1, Achi) + kron(H2, Aabs));
        rhs = kron(Q1, Aplus*rhoL) + kron(Q2, -Aminus*rhoR);

    case {'central', 'central-local', 'local-central'}
        AcentralRho = localCentralRelativeOperatorRho(p, size(Phi, 1));
        [AplusBoundary, AminusBoundary] = splitByCharacteristicSign(AcentralRho);
        AabsBoundary = AplusBoundary - AminusBoundary;
        A = p.Q_diff*(kron(H1, AcentralRho) + kron(H2Boundary, AabsBoundary));
        rhs = kron(Q1, AplusBoundary*rhoL) + kron(Q2, -AminusBoundary*rhoR);

    case {'rusanov', 'lax-friedrichs', 'lf', 'rusanov-local', 'local-rusanov'}
        % Sparse local rho-basis transport operator
        AcentralRho = localCentralRelativeOperatorRho(p, size(Phi, 1));
        nRho = size(AcentralRho, 1);

        % Rusanov / Lax-Friedrichs strength
        if isfield(mat.dg.params, 'thetaLF')
            thetaLF = mat.dg.params.thetaLF;
        elseif isfield(mat.dg.params, 'rho_theta_lf')
            thetaLF = mat.dg.params.rho_theta_lf;
        else
            thetaLF = 0.1;  % safe default for testing
        end

        if isfield(mat.dg.params, 'alphaLF')
            alphaLF = mat.dg.params.alphaLF;
        elseif isfield(mat.dg.params, 'rho_alpha_lf')
            alphaLF = mat.dg.params.rho_alpha_lf;
        else
            alphaLF = max(abs(real(lambda)));
            if alphaLF == 0
                alphaLF = max(abs(lambda));
            end
        end

        betaLF = thetaLF * alphaLF;
        % Sparse scalar LF approximation:
        % |A| ≈ alphaLF * I
        AabsLF = betaLF * speye(nRho);
        % Characteristic upwind only at physical chi-boundaries.
        % This is dense in general, but only appears in boundary blocks.
        [AplusBoundary, AminusBoundary] = splitByCharacteristicSign(AcentralRho);
        AabsBoundary = AplusBoundary - AminusBoundary;
        % H2 contains the absolute/upwind part for all faces.
        % H2Boundary contains only the two physical boundary faces.
        % Therefore H2Interior applies LF stabilization only to interior faces.
        H2Interior = H2 - H2Boundary;
        A = p.Q_diff*( ...
              kron(H1, AcentralRho) ...
            + kron(H2Interior, AabsLF) ...
            + kron(H2Boundary, AabsBoundary) );

        rhs = kron(Q1, AplusBoundary*rhoL) ...
            + kron(Q2, -AminusBoundary*rhoR);

    case {'rusanov-boundary-upwind'}
        % Sparse local rho-basis transport operator
        AcentralRho = localCentralRelativeOperatorRho(p, size(Phi, 1));
        nRho = size(AcentralRho, 1);
    
        % ------------------------------------------------------------
        % Rusanov / Lax-Friedrichs strength
        % ------------------------------------------------------------
        if isfield(mat.dg.params, 'thetaLF')
            thetaLF = mat.dg.params.thetaLF;
        elseif isfield(mat.dg.params, 'rho_theta_lf')
            thetaLF = mat.dg.params.rho_theta_lf;
        else
            thetaLF = 0.01;
        end
    
        if isfield(mat.dg.params, 'alphaLF')
            alphaLF = mat.dg.params.alphaLF;
        elseif isfield(mat.dg.params, 'rho_alpha_lf')
            alphaLF = mat.dg.params.rho_alpha_lf;
        else
            alphaLF = max(abs(real(lambda)));
            if alphaLF == 0
                alphaLF = max(abs(lambda));
            end
        end
    
        betaLF = thetaLF * alphaLF;
    
        % Sparse scalar LF approximation:
        % |A| ≈ betaLF * I
        AabsLF = betaLF * speye(nRho);
    
        % ------------------------------------------------------------
        % Characteristic upwind matrix for the boundary layer
        % ------------------------------------------------------------
        [AplusBoundary, AminusBoundary] = splitByCharacteristicSign(AcentralRho);
        AabsBoundary = AplusBoundary - AminusBoundary;
    
        % Optional: store as sparse. It may still be numerically dense,
        % but only in the selected boundary-layer blocks.
        AabsBoundary = sparse(AabsBoundary);
        AplusBoundary = sparse(AplusBoundary);
        AminusBoundary = sparse(AminusBoundary);
    
        % ------------------------------------------------------------
        % Number of chi-cells on each side where characteristic upwind
        % should be used.
        %
        % nUpwindLayer = 0: only physical boundary faces, as before
        % nUpwindLayer = 2: first/last two chi-cell equations use
        %                   characteristic upwind in the absolute flux part
        % ------------------------------------------------------------
        if isfield(mat.dg.params, 'rho_upwind_layers')
            nUpwindLayer = mat.dg.params.rho_upwind_layers;
        elseif isfield(mat.dg.params, 'upwindBoundaryLayers')
            nUpwindLayer = mat.dg.params.upwindBoundaryLayers;
        else
            nUpwindLayer = 0;
        end
    
        if nUpwindLayer <= 0
            % Previous behavior: characteristic upwind only on physical
            % left/right boundary faces.
            H2Upwind = H2Boundary;
        else
            % New behavior: characteristic upwind in the first/last
            % nUpwindLayer chi-cell equations.
            H2Upwind = boundaryLayerRows(H2, nUpwindLayer);
        end
    
        % Remaining interior part gets sparse Rusanov/LF stabilization
        H2Rusanov = H2 - H2Upwind;
    
        A = p.Q_diff*( ...
              kron(H1, AcentralRho) ...
            + kron(H2Rusanov, AabsLF) ...
            + kron(H2Upwind, AabsBoundary) );
    
        % RHS still only comes from physical inflow boundaries
        rhs = kron(Q1, AplusBoundary*rhoL) ...
            + kron(Q2, -AminusBoundary*rhoR);
    case {'central-projected', 'projected-central'}
        AcentralRho = Achi;
        A = p.Q_diff*(kron(H1, AcentralRho) + kron(H2Boundary, Aabs));
        rhs = kron(Q1, Aplus*rhoL) + kron(Q2, -Aminus*rhoR);
    otherwise
        error('Unknown rho flux "%s". Use "central" or "matrix-upwind".', fluxType);
end

rhs = p.Q_diff*rhs;

%% Flux information for afterwards control
fluxInfo = struct;
fluxInfo.type = fluxType;
fluxInfo.lambdaRange = [min(real(lambda)), max(real(lambda))];
fluxInfo.phiGramDiagRange = [min(real(diag(Phi'*Phi))), max(real(diag(Phi'*Phi)))];
if exist('PhiLeft', 'var')
    fluxInfo.phiLeftResidual = norm(PhiLeft*Phi - eye(size(Phi, 2)), 'fro');
end
if exist('Aplus', 'var')
    fluxInfo.Aplus = Aplus;
    fluxInfo.Aminus = Aminus;
end
if any(strcmpi(fluxType, {'central', 'central-local', 'local-central', 'central-projected', 'projected-central'}))
    fluxInfo.A = AcentralRho;
    fluxInfo.Acentral = AcentralRho;
else
    fluxInfo.A = Achi;
end
fluxInfo.rhoL = rhoL;
fluxInfo.rhoR = rhoR;
end
%% 



%% Helper Functions
function H2Boundary = boundaryUpwindChiMatrix(p, detJx, K2L, K2R)
alphaL = zeros(p.N_chi, 1);
alphaR = zeros(p.N_chi, 1);
alphaL(1) = 1;
alphaR(end) = 1;

H2Boundary = kron(spdiags((1./detJx).*alphaL, 0, p.N_chi, p.N_chi), K2L) ...
           + kron(spdiags((1./detJx).*alphaR, 0, p.N_chi, p.N_chi), K2R);
end

function PhiLeft = leftInversePhi(Phi)
gram = Phi'*Phi;
PhiLeft = gram\Phi';
end


function [Aplus, Aminus] = splitByCharacteristicSign(A)
[V, Lambda] = eig(full(A));
lambda = diag(Lambda);
lambdaPlus = (real(lambda) > 0).*lambda;
lambdaMinus = (real(lambda) < 0).*lambda;

if norm(A -A', 'fro') <= 1e-12*max(1, norm(A, 'fro'))
    Aplus = V*diag(lambdaPlus)*V';
    Aminus = V*diag(lambdaMinus)*V';
else
    Vinv = V\eye(size(V));
    Aplus = V*diag(lambdaPlus)*Vinv;
    Aminus = V*diag(lambdaMinus)*Vinv;
end
end

% function rhs = centralBoundaryRhs(Q1, Q2, AcentralRho, rhoL, rhoR)
% rhs = kron(Q1, 0.5*AcentralRho*rhoL) - kron(Q2, 0.5*AcentralRho*rhoR);
% end

function AcentralRho = localCentralRelativeOperatorRho(p, nRho)
% Use the local relative-coordinate k-operator for central fluxes. The
% rho basis stores relative-coordinate values, so k maps to -1i*d/dxi. The
% projected central operator is generally denser and is kept in the separate
% 'central-projected' diagnostic mode.
if isfield(p, 'Dy') && ~isscalar(p.Dy) && isequal(size(p.Dy), [nRho, nRho])
    AcentralRho = -1i*sparse(p.Dy);
elseif p.doDG
    AcentralRho = -1i*sparse(get_Discretization(p)/p.delta_xi/2);
elseif isfield(p, 'Dy') && isequal(size(p.Dy), [nRho, nRho])
    AcentralRho = -1i*sparse(p.Dy);
else
    error('No sparse rho-basis central xi-operator available for size %d.', nRho);
end
end

function Hlayer = boundaryLayerRows(H, nLayer)
%BOUNDARYLAYERROWS Select first/last nLayer chi-cell equation rows.
%
% This creates a sparse matrix with the same entries as H, but only
% in the first and last nLayer block rows in chi. It is used to apply
% characteristic upwind only near the physical boundaries, while the
% interior keeps the sparse Rusanov/LF stabilization.

N = size(H, 1);

nLayer = max(0, round(nLayer));
nLayer = min(nLayer, floor(N/2));

if nLayer == 0
    Hlayer = spalloc(N, N, 0);
    return;
end

rowMask = false(N, 1);
rowMask(1:nLayer) = true;
rowMask(N-nLayer+1:N) = true;

R = spdiags(double(rowMask), 0, N, N);

Hlayer = R * H;
Hlayer = sparse(Hlayer);
end
