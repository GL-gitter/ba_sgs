function DG = solve_transport_DG_2D_transient_SC_new_rho(mat, INIT)
%SOLVE_TRANSPORT_DG_2D_TRANSIENT_SC_NEW_RHO RK4 transient rho-basis variant.
%
% INIT is expected to come from solve_transport_DG_2D_rho with
% mat.dg.params.st_init = true, so INIT.rho stores rho-basis vectors.

addpath DG_code/core/

h = 6.626E-34;
m0 = 9.109E-31;
hb = h/2/pi;

N_K_chi = mat.dg.params.N_K_chi;
N_chi = ((mat.Nx-1)/(N_K_chi-1));
mat.Np = N_chi*N_K_chi;

deg_factor = mat.deg_factor;
n_of_modes = mat.n_of_modes;

Vsub_INIT = INIT.Vsub;
rhs_INIT = INIT.rhs;
rho_INIT = INIT.rho;
p = INIT.p;
Vbi = INIT.V_total;

N_rho = size(p.phi, 1);
rho_pre = zeros(n_of_modes, p.N_chi*p.N_K_chi, N_rho);
for IM = 1:n_of_modes
    rho_pre(IM, :, :) = reshape(squeeze(rho_INIT{IM}), N_rho, p.N_chi*p.N_K_chi).';
end

Vsub_mod = zeros(n_of_modes, p.N_chi*p.N_K_chi, mat.Ny);
Vsub_mod_pre = zeros(n_of_modes, p.N_chi*p.N_K_chi, mat.Ny);
m_mod = zeros(n_of_modes, p.N_chi*p.N_K_chi, 1);
n_pre = zeros(n_of_modes, p.N_chi*p.N_K_chi, mat.Ny);

for IM = 1:n_of_modes
    m = squeeze(mat.me_x(1, :, ceil(mat.Ny/2)))*m0;

    for Ichi = 1:p.N_chi
        chiRange = (Ichi-1)*p.N_K_chi+1:Ichi*p.N_K_chi;
        srcRange = (Ichi-1)*(p.N_K_chi-1)+1:(Ichi)*(p.N_K_chi-1)+1;

        V_temp = reshape(Vsub_INIT(1, srcRange, IM, :), p.N_K_chi, mat.Ny);
        Vsub_mod(IM, chiRange, :) = V_temp;
        Vsub_mod_pre(IM, chiRange, :) = V_temp;

        m_temp = reshape(m(srcRange), p.N_K_chi, 1);
        m_mod(IM, chiRange) = m_temp;
    end

    if mod(N_rho, 2) == 0
        n_center = squeeze(real(squeeze(rho_pre(IM, :, N_rho/2+1)) + squeeze(rho_pre(IM, :, N_rho/2)))/2);
    else
        n_center = squeeze(real(rho_pre(IM, :, (N_rho-1)/2+1)));
    end
    n_pre(IM, :, :) = kron(ones(1, mat.Ny), n_center.').*reshape(abs(Vsub_mod_pre(IM, :, :)).^2, p.N_chi*p.N_K_chi, mat.Ny)/mat.dy/1E-9;
end

n_t_pre = squeeze(sum(n_pre, 1));
m_mod = squeeze(m_mod(1, :, :)).';

n_over_t = zeros(p.N_chi*p.N_K_chi, mat.Ny, mat.Nt);
j_over_t = zeros(p.N_chi*p.N_K_chi, mat.Nt);
n_over_t_mode = zeros(n_of_modes, p.N_chi*p.N_K_chi, mat.Ny, mat.Nt);
j_over_t_mode = zeros(n_of_modes, p.N_chi*p.N_K_chi, mat.Nt);
fid = figure('name', "Dichte rho");

Nt = mat.Nt;
dt = mat.dt;

n_sp = zeros(mat.Nx, mat.Ny);
MapM = [0.5, 0, 0; 0, 1, 0; 0, 0, 0.5];
Mapping = zeros(mat.Nx, mat.Np);
for IN = 1:N_chi
    Mapping((IN-1)*(N_K_chi-1)+1:IN*(N_K_chi-1)+1, (IN-1)*N_K_chi+1:IN*N_K_chi) = MapM;
end
Mapping(1, 1) = 1;
Mapping(mat.Nx, mat.Np) = 1;

IV = 1;
A = cell(1, n_of_modes);

for tstep = 1:Nt
    if mod((tstep-1), 50) == 0
        display(['Zeitschritt ', num2str(tstep), '/', num2str(Nt)])
    end

    if (mod((tstep-1), 20) == 0 || (tstep < 50))
        Vref = mat.V;
        mat.V = Vref + Vbi;
        [EM, VM] = solve_subbands(mat); %#ok<ASGLU>
        mat.V = Vref;
        [Ef0_L, Ef0_R] = solve_fermi(EM, mat);
    end

    for IM = 1:n_of_modes
        if (mod((tstep-1), 20) == 0 || (tstep < 50))
            m = squeeze(mat.me_x(IV, :, ceil(mat.Ny/2)))*m0;
            p = initParams(mat, Ef0_L-mat.Vs, Ef0_R-mat.Vd, squeeze(EM(IV, :, IM)), m, deg_factor(IV));
            [A{IM}, ~] = get_SysM_rho(mat, p, squeeze(EM(IV, :, IM)));
            A{IM} = mat.dg.params.MPC*A{IM};
        end

        k1 = rhs_INIT{IM} - A{IM}*rho_INIT{IM};
        k2 = rhs_INIT{IM} - A{IM}*(rho_INIT{IM} + dt/2*k1);
        k3 = rhs_INIT{IM} - A{IM}*(rho_INIT{IM} + dt/2*k2);
        k4 = rhs_INIT{IM} - A{IM}*(rho_INIT{IM} + dt*k3);

        rho_INIT{IM} = rho_INIT{IM} + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
        result_matrix = reshape(squeeze(rho_INIT{IM}), N_rho, p.N_chi*p.N_K_chi).';

        if mod(N_rho, 2) == 0
            n_over_t_mode(IM, :, :, tstep) = kron(ones(1, mat.Ny), real((result_matrix(:, N_rho/2+1) + result_matrix(:, N_rho/2))/2)).*reshape(abs(Vsub_mod(IM, :, :)).^2, p.N_chi*p.N_K_chi, mat.Ny)/mat.dy/1E-9;
            j_over_t_mode(IM, :, tstep) = p.q*hb./m_mod.*imag(result_matrix(:, N_rho/2+1) - result_matrix(:, N_rho/2))/p.delta_xi;
        else
            n_over_t_mode(IM, :, :, tstep) = kron(ones(1, mat.Ny), real(result_matrix(:, (N_rho-1)/2+1))).*reshape(abs(Vsub_mod(IM, :, :)).^2, p.N_chi*p.N_K_chi, mat.Ny)/mat.dy/1E-9;
            j_over_t_mode(IM, :, tstep) = p.q*hb./m_mod.*imag(result_matrix(:, (N_rho-1)/2+2) - result_matrix(:, (N_rho-1)/2))/p.delta_xi;
        end
    end

    if (n_of_modes > 1)
        n_over_t(:, :, tstep) = sum(squeeze(n_over_t_mode(:, :, :, tstep)), 1);
        j_over_t(:, tstep) = sum(squeeze(j_over_t_mode(:, :, tstep)), 1);
    else
        n_over_t(:, :, tstep) = squeeze(n_over_t_mode(:, :, :, tstep));
        j_over_t(:, tstep) = squeeze(j_over_t_mode(:, :, tstep));
    end

    if mod((tstep-1), 50) == 0
        figure(fid)
        title(['t=', num2str(mat.t(tstep)), 'fs']);
        plot(p.chi_Koord, real(squeeze(n_t_pre(:, ceil(mat.Ny/2)))), 'r');
        hold on
        plot(p.chi_Koord, real(squeeze(n_over_t(:, ceil(mat.Ny/2), tstep))), 'b');
        hold off
        drawnow limitrate
    end

    if (mod(tstep, 20) == 0 || (tstep < 50))
        for IY = 1:mat.Ny
            n_sp(:, IY) = Mapping*squeeze(n_over_t(:, IY, tstep));
        end
        Vbi_IN = INIT.V_total;
        Vbi = solve_Poisson2D_transient(mat, Vbi_IN, INIT.Ef0, tstep, n_sp);
    end
end

DG.basis = 'rho';
DG.n_over_t = n_over_t(:, :, 1:50:end);
DG.j_over_t = j_over_t(:, 1:5:end);
DG.rho_INIT = rho_INIT;
DG.chi = p.chi_Koord;
DG.xi = p.xi_Rechengebiet_Elem;
end
