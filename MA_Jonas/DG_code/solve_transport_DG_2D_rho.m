function DG = solve_transport_DG_2D_rho(mat, Esub, Vsub, EfL, EfR)
%SOLVE_TRANSPORT_DG_2D_RHO Stationary 1D Wigner-DG solve in rho basis.
%
% This is a parallel implementation to solve_transport_DG_2D. The existing
% phase/eigen-basis solver is not modified.

h = 6.626E-34;
m0 = 9.109E-31;
hb = h/2/pi;

addpath DG_code/core/

DG = struct;

n_of_valleys = mat.n_of_valleys;
n_of_modes = mat.n_of_modes;
deg_factor = mat.deg_factor;

N_K_chi = mat.dg.params.N_K_chi;
N_chi = ((mat.Nx-1)/(N_K_chi-1));
N_chi_dof = N_chi*N_K_chi;

p0 = initParams(mat, EfL, EfR, squeeze(Esub(1, :, 1)), squeeze(mat.me_x(1, :, ceil(mat.Ny/2)))*m0, deg_factor(1));
N_rho = size(p0.phi, 1);
RHO_DG = zeros(n_of_valleys, n_of_modes, N_chi_dof, N_rho);

n_mode = zeros(n_of_valleys, n_of_modes, N_chi_dof, mat.Ny);
j_mode = zeros(n_of_valleys, n_of_modes, N_chi_dof, 1);
IN_max = size(n_mode, 3);

p_top = repmat(p0, n_of_valleys, n_of_modes);
for IV = 1:n_of_valleys
    for IM = 1:n_of_modes
        p = initParams(mat, EfL, EfR, squeeze(Esub(IV, :, IM)), squeeze(mat.me_x(IV, :, ceil(mat.Ny/2)))*m0, deg_factor(IV));
        p_top(IV, IM) = p;
    end
end

if (mat.dg.params.st_init == true)
    Atemp = cell(1, n_of_modes);
    rhstemp = cell(1, n_of_modes);
    rhotemp = cell(1, n_of_modes);
end
info = cell(n_of_valleys, n_of_modes);

for IV = 1:n_of_valleys
    for IM = 1:n_of_modes
        m = squeeze(mat.me_x(IV, :, ceil(mat.Ny/2)))*m0;
        p = p_top(IV, IM);

        Vsub_mod = zeros(p.N_chi*p.N_K_chi, mat.Ny);
        m_mod = zeros(p.N_chi*p.N_K_chi, 1);

        for Nchi = 1:p.N_chi
            V_temp = reshape(Vsub(IV, (Nchi-1)*(p.N_K_chi-1)+1:(Nchi)*(p.N_K_chi-1)+1, IM, :), p.N_K_chi, mat.Ny);
            Vsub_mod((Nchi-1)*p.N_K_chi+1:Nchi*p.N_K_chi, :) = V_temp;

            m_temp = reshape(m((Nchi-1)*(p.N_K_chi-1)+1:(Nchi)*(p.N_K_chi-1)+1), p.N_K_chi, 1);
            m_mod((Nchi-1)*p.N_K_chi+1:Nchi*p.N_K_chi) = m_temp;
        end
        m_mod = m_mod.';

        [A, rhs, sysInfo] = get_SysM_rho(mat, p, squeeze(Esub(IV, :, IM)));

        if (mat.dg.params.MPC ~= 1)
            A = mat.dg.params.MPC*A;
            rhs = mat.dg.params.MPC*rhs;
        end

        rho = A\rhs;

        if (mat.dg.params.st_init == true)
            Atemp{IM} = A;
            rhstemp{IM} = rhs;
            rhotemp{IM} = rho;
        end
        info{IV, IM} = sysInfo;

        rho_chi_xi = reshape(rho, N_rho, p.N_chi*p.N_K_chi).';
        rho_xi_chi = rho_chi_xi.';
        X = rho_xi_chi;

        if (mat.dg.params.sl == 1)
            legendre_limited = zeros(N_rho, p.N_chi*p.N_K_chi);
            for IN = 1:p.N_chi*p.N_K_chi
                legendre_limited(:, IN) = apply_legendre_limiter(real(rho_xi_chi(:, IN)), p.N_K_xi, 1);
            end
            X = legendre_limited;
        end

        if (p.doFV == true)
            if mod(p.N_xi, 2) == 1
                n_DG = real(X((p.N_xi-1)/2+1, :))';
                j_mode(IV, IM, :) = (p.q*hb./m_mod.*imag(rho_xi_chi((p.N_xi-1)/2+2, :) - rho_xi_chi((p.N_xi-1)/2, :))/2/p.delta_xi)';
            else
                n_DG = ((real(X(p.N_xi/2, :)) + real(X(p.N_xi/2+1, :)))/2)';
                j_mode(IV, IM, :) = (p.q*hb./m_mod.*imag(rho_xi_chi(p.N_xi/2+1, :) - rho_xi_chi(p.N_xi/2, :))/p.delta_xi)';
            end
        elseif (p.doDG == true)
            N_xi_total = p.N_xi*p.N_K_xi;
            if mod(N_xi_total, 2) == 1
                n_DG = real(X((N_xi_total-1)/2+1, :))';
                j_mode(IV, IM, :) = (p.q*hb./m_mod.*imag(rho_xi_chi((N_xi_total-1)/2-2, :) - rho_xi_chi((N_xi_total-1)/2, :))/2/p.delta_xi)';
            else
                n_DG = ((real(X(N_xi_total/2, :)) + real(X(N_xi_total/2+1, :)))/2)';
                j_mode(IV, IM, :) = (p.q*hb./m_mod.*imag(rho_xi_chi(N_xi_total/2-1, :) - rho_xi_chi(N_xi_total/2+2, :))/p.delta_xi)';
            end
        else
            error('solve_transport_DG_2D_rho currently supports FV and DG in xi.');
        end

        for IN = 1:IN_max
            n_mode(IV, IM, IN, :, :) = (squeeze(n_DG(IN, 1)).*squeeze(abs(Vsub_mod(IN, :))).^2/mat.dy/1E-9);
        end

        RHO_DG(IV, IM, :, :) = rho_chi_xi;
    end
end

nT = squeeze(sum(n_mode, [1 2]));
jT = squeeze(sum(j_mode, [1 2]));

if (mat.dg.params.st_init == true)
    DG.A = Atemp;
    DG.rhs = rhstemp;
    DG.rho = rhotemp;
end

DG.basis = 'rho';
DG.flux = 'matrix-upwind';
if isfield(mat.dg.params, 'rho_flux')
    DG.flux = mat.dg.params.rho_flux;
end
DG.info = info;
DG.p = p_top(1, 1);
DG.n = nT;
DG.j = jT;
DG.RHO = RHO_DG;
DG.chi = p_top(1, 1).chi_Koord;
DG.xi = p_top(1, 1).xi_Rechengebiet_Elem;
end

