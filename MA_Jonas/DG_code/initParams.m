function p = initParams(mat, EfL, EfR, V, m, deg_factor)

p = struct;

p.permute_shell = mat.dg.params.permute_shell;

p.deg_factor_temp = deg_factor;
p.doMSA_3D = false;
p.doMSA_2D = false;

%% Choose MSA 3D, MSA 2D or 1D
if (strcmpi(mat.type,'GAAFET')==true)
    p.doMSA_3D = true;
elseif ((strcmpi(mat.type,'DGFET')==true) || (strcmpi(mat.type,'DGFET-type-I')==true))
    p.doMSA_2D = true;
end

%% Chose Dimension
p.do1D = true;         % DG algoirthm with 1D finite elements
p.do2D = false;          % DG algorithm with 2D finite elements

%% Choose Simulation Case
p.doSteadyState = true;
p.doTransient = false;
p.doIV = false;
p.doPoisson = false;

if (p.doMSA_3D == true)
    p.doPoisson   = false;
    p.do1D        = true;
    p.doIV        = false;
    p.doTransient = false;
end

%% Choose numerical algorithm for xi
if (strcmpi(mat.dg.params.Xi_solver, 'FV')==true)
    p.doDG = false;       
    p.doFEM = false;
    p.doFV = true;
elseif (strcmpi(mat.dg.params.Xi_solver, 'DG')==true)
    p.doDG = true;       
    p.doFEM = false;
    p.doFV = false;
end

%% Set best fermi distribution acc. to simulation case
if (p.doPoisson == true && p.doMSA == false)
    p.choose_fermi = 'LS';
elseif (p.doMSA_3D == true)
    p.choose_fermi = 'DG_MSA_3D';
elseif (p.doMSA_2D == true)
    p.choose_fermi = 'DG_MSA_2D';
else
    p.choose_fermi = 'A';
end

%% Numerical Constants
p.N_K_chi = mat.dg.params.N_K_chi;
p.N_chi   = ((mat.Nx-1)/(p.N_K_chi-1));    %
p.N_xi    = mat.dg.params.N_xi;       % even for FV in xi and odd for FEM in xi 
p.Nqx     = 10;                                          % Ordnung der Gauß-Lobatto-Quadratur (entspr. Anz. Stützstellen in der Einheitszelle und Anz. Gewichte)
p.Nk     = mat.dg.params.N;
p.BT     = mat.dg.params.Q;

if (p.doFV == true)
    p.N_K_xi = 1;           % Do not change this value!!!
    p.Ny = p.N_xi+1;
    p.dim = p.N_chi*p.N_xi*p.N_K_chi;
    p.Nk_K = 1;
elseif (p.doDG == true)
    p.N_K_xi = mat.dg.params.N_K_xi;           % Set number of nodes for Xi-domain
    p.Ny = p.Nqx*p.N_xi;
    p.dim = p.N_chi*p.N_K_chi*p.N_xi*p.N_K_xi;
    p.Nk_K   = mat.dg.params.Nk_K;
elseif (p.doFEM == true)
    p.N_K_xi = 2;           % Set number of nodes for Xi-domain
    p.Ny = p.Nqx*(p.N_xi*(p.N_K_xi-1)+1);
    p.dim = p.N_chi*p.N_K_chi*(p.N_xi*(p.N_K_xi-1)+1);
end

%% Variables concerning IV
if (p.doIV == true)
    p.U0 = 0.0;
    p.U = 0.35;
    p.NU = 5;
    p.U_IV = linspace(p.U0,p.U,p.NU);
end

%% Variables concerning transient simulation
if (p.doTransient == true)
    p.U = 0.2;
    p.t = 0;                                                          % Startzeit
    p.TFinal = 1e-13;                                                  % Endzeit
    p.dt = 1e-17;
    p.tsteps = ceil(p.TFinal/p.dt);
    p.time = linspace(p.t, p.TFinal, p.tsteps/10);
    p.dt = p.TFinal/p.tsteps;
end

%% Physical Constants
p.hquer = 1.054571817e-34;                        % Reduziertes Plancksches Wirkungsquantum in [Js]
p.m_elektron = 0.041*9.1094e-31;                         % Masse eines Elektrons in [kg]
p.q = 1.602176634e-19;                            % Ladung eines Elektrons in [C]
p.kB = 1.38064852e-23;
p.T = 300;
p.Q_diff = p.hquer/p.m_elektron;                              % Einheit [Js/kg]
p.Q_drift = p.q/p.hquer;                                      % Einheit [C/Js]
p.U0 = 0.0;                                            %Startspannung für rampTime
p.U = 0.0;                          % Spannung für thermisches Nichtgleichgewicht
p.L_xi = mat.dg.params.Ly;        %Länge in xi-Richtung?
p.L_chi = mat.L_x*1e-9;           %60e-9; 112.5e-9        %Länge in chi-Richtung?
p.L_d = 24e-9;            %24e-9;  51.5e-9;      %Länge des nicht dotierten Bereichs?
p.bar_g = 3e-9;          %5e-9;  3e-9;       % entspr. L_b aus Abb. 7.1
p.bar_w = 3.5e-9;          %5e-9;  3.5e-9;       % entspr. L_w aus Abb. 7.1
p.alpha = 0;
% p.rampTime = 100e-15;
p.eps_0 = 8.854e-12;
p.eps_r_GaAs = 12.9;
p.eps_r_AlGaAs = 11.9628;
p.beta = 1/(p.kB*p.T);

p.N = p.N_K_chi+1;
p.Np = (p.N-1)*p.N_chi+1;
p.xg = linspace(0, p.L_chi, p.Np)';
p.dx = (p.xg(2)-p.xg(1));
p.Nd_V = ones(1, p.Np)*1e24;      %1e24;
p.Nd_V(p.xg < p.L_chi/2+p.L_d/2 & p.xg > p.L_chi/2-p.L_d/2) = 0;
if (p.doPoisson == false)
    p.eps_V = ones(1, p.N_chi*p.Nqx)*p.eps_r_GaAs*p.eps_0;                      %Warum Nqx hier?
    p.eps_V(p.xg < p.L_chi/2+p.bar_g+p.bar_w/2) = p.eps_r_AlGaAs*p.eps_0;
    p.eps_V(p.xg < p.L_chi/2-p.bar_g-p.bar_w/2) = p.eps_r_GaAs*p.eps_0;
    p.eps_V(abs(p.xg) >= p.L_chi/2-p.bar_w/2 & abs(p.xg) <= p.L_chi/2+p.bar_w/2) = p.eps_r_GaAs*p.eps_0;
    p.eps_V_he = [p.eps_V(1), p.eps_V, p.eps_V(p.N_chi*p.Nqx)];
    p.eps_V = (p.eps_V_he(1:end-1)+p.eps_V_he(2:end))/2;
else
    p.eps_V = ones(1, p.Np)*p.eps_r_GaAs*p.eps_0;
    p.eps_V(p.xg < p.L_chi/2+p.bar_g+p.bar_w/2) = p.eps_r_AlGaAs*p.eps_0;
    p.eps_V(p.xg < p.L_chi/2-p.bar_g-p.bar_w/2) = p.eps_r_GaAs*p.eps_0;
    p.eps_V(abs(p.xg) >= p.L_chi/2-p.bar_w/2 & abs(p.xg) <= p.L_chi/2+p.bar_w/2) = p.eps_r_GaAs*p.eps_0;
    p.eps_V_he = [p.eps_V(1), p.eps_V, p.eps_V(p.Np)];
    p.eps_V = (p.eps_V_he(1:end-1)+p.eps_V_he(2:end))/2;
end
p.err = 1;
p.V_ref = 1/(p.q*p.beta);

if (p.doPoisson == true && p.doIV == false || p.doTransient == true || p.doSteadyState == true)
    p.U1 = zeros(size(p.xg));
    p.U1(p.xg < p.L_chi/2+p.L_d/2 & p.xg > p.L_chi/2-p.L_d/2) = linspace(0, p.U, sum(p.xg < p.L_chi/2+p.L_d/2 & p.xg > p.L_chi/2-p.L_d/2));
    p.U1(p.xg >= p.L_chi/2+p.L_d/2) = p.U;
    p.UU = p.U1;
end

if (p.do1D == true)
    %% Discretize chi-domain
    p.delta_chi = p.L_chi/p.N_chi;
    p.disk_Rechengebiet1 = linspace(0, p.L_chi, p.N_chi+1);
    for a = 1 : p.N_chi
        p.chi_Koord(1, p.N_K_chi*(a-1)+1:p.N_K_chi*a) = linspace(p.disk_Rechengebiet1(a), p.disk_Rechengebiet1(a+1), p.N_K_chi);
    end

    %% Get Potential
    p.V_pot = get_V(p.chi_Koord,p)';

    %% Discretize xi-domain
    if (p.doDG == true)
        p.xi_Rechengebiet = linspace(-p.L_xi/2, p.L_xi/2, p.N_xi+1);
        for a = 1 : p.N_xi
            p.xi_Koord(1, p.N_K_xi*(a-1)+1:p.N_K_xi*a) = linspace(p.xi_Rechengebiet(a), p.xi_Rechengebiet(a+1), p.N_K_xi);
        end
        p.delta_xi = diff(p.xi_Rechengebiet);
        p.delta_xi = p.delta_xi(1);
        p.xi_Rechengebiet_Elem = linspace(-(p.L_xi-p.delta_xi)/2, (p.L_xi-p.delta_xi)/2, p.N_xi*p.N_K_xi);
    elseif (p.doFEM == true)
        p.xi_Rechengebiet = linspace(-p.L_xi/2, p.L_xi/2, p.N_xi*(p.N_K_xi-1)+1);
        p.delta_xi = diff(p.xi_Rechengebiet);
        p.delta_xi = p.delta_xi(1);
        p.xi_Rechengebiet_Elem = linspace(-(p.L_xi-p.delta_xi)/2, (p.L_xi-p.delta_xi)/2, p.N_xi+1);
    elseif (p.doFV == true)
        p.xi_Rechengebiet = linspace(-p.L_xi/2, p.L_xi/2, p.N_xi);
        p.delta_xi = diff(p.xi_Rechengebiet);
        p.delta_xi = p.delta_xi(1);
        p.xi_Rechengebiet_Elem = linspace(-(p.L_xi-p.delta_xi)/2, (p.L_xi-p.delta_xi)/2, p.N_xi);
    end

    %% Get DG-Mass matrix in chi
    [p.x, p.w] = JacobiGL(0, 0, p.Nqx-1);
    p.r = JacobiGL(0, 0, p.N_K_chi-1);
    p.V_tilde = Vandermonde1D(p.N_K_chi-1, p.x);
    p.V = Vandermonde1D(p.N_K_chi-1, p.r);
    p.invV = inv(p.V);
    p.M = p.invV.'*p.invV;

    %% Get DG/FEM-Mass matrix in xi
    [p.x_xi, p.w_xi] = JacobiGL(0, 0, p.Nqx-1);
    p.r_xi = JacobiGL(0, 0, p.N_K_xi-1);
    p.V_tilde_xi = Vandermonde1D(p.N_K_xi-1, p.x);
    p.V_xi = Vandermonde1D(p.N_K_xi-1, p.r_xi);
    p.M_xi = inv(p.V_xi).'*inv(p.V_xi);

    if (p.doDG == true)
        p.M_xiCell = repmat({p.M_xi},1,p.N_xi);
        p.M_xi_glob = blkdiag(p.M_xiCell{:});
    elseif (p.doFEM == true)
        p.M_xi_glob = zeros(p.N_xi*(p.N_K_xi-1)+1,p.N_xi*(p.N_K_xi-1)+1);
        p.M_xi_glob(1:p.N_K_xi,1:p.N_K_xi) = p.M_xi;
        for aa = 1:p.N_xi-1
            p.M_xi_glob(aa+1:aa+p.N_K_xi, aa+1:aa+p.N_K_xi)=[p.M_xi_glob(aa+1,aa+1)+p.M_xi(1,1), p.M_xi(1,p.N_K_xi);
                p.M_xi(p.N_K_xi,1),                 p.M_xi(p.N_K_xi,p.N_K_xi)];
        end
    end

    %% Get Eigenvectors and Eigenvalues
    [p.Dy, p.phi,p.D,p.k,p.dk,p.Q] = Basisfunktionen_Q_m(p, p.N_xi*p.N_K_xi);      % Eigenvektoren und Eigenwerte
    % p.Wellenintervall = diag(p.D);

    %% Get Fermi distribution
    if (p.doFEM == true)
        [p.verteilung] = verteilung(p, p.N_xi*(p.N_K_xi-1)+1, p.phi);
    else
        [p.verteilung_l, p.verteilung_r] = verteilung(mat, p, p.N_xi*p.N_K_xi, p.phi, EfL, EfR, V, m);                       % Gaussche Verteilungsfunktion
    end









    

    %% 2D


elseif (p.do2D == true)
    % 2D numerical constants
    p.N_xi2D = 30;
    p.N_chi2D = 20;
    p.N2Dp = 2;
    p.N2D = p.N_xi2D*p.N_chi2D;
    p.K = zeros(p.N2Dp*p.N_xi2D,p.N2Dp*p.N_chi2D);
    for i = 1:p.N_xi2D
        for j = 1:p.N_chi2D
            p.K((i-1)*p.N2Dp+1,(j-1)*p.N2Dp+1:(j-1)*p.N2Dp+2) = (i-1)*p.N_chi2D*p.N2Dp*p.N2Dp+(j-1)*4+1:(i-1)*p.N_chi2D*p.N2Dp*p.N2Dp+(j-1)*4+2;
            p.K((i-1)*p.N2Dp+2,(j-1)*p.N2Dp+1:(j-1)*p.N2Dp+2) = fliplr((i-1)*p.N_chi2D*p.N2Dp*p.N2Dp+(j-1)*4+3:(i-1)*p.N_chi2D*p.N2Dp*p.N2Dp+(j-1)*4+4);
        end
    end
    p.Ky = zeros(p.N2Dp*p.N_xi2D,p.N2Dp);
    for i = 1:p.N_xi2D
        p.Ky((i-1)*p.N2Dp+1,1:2) = (i-1)*p.N2Dp*p.N2Dp+1:(i-1)*p.N2Dp*p.N2Dp+2;
        p.Ky((i-1)*p.N2Dp+2,1:2) = fliplr((i-1)*p.N2Dp*p.N2Dp+3:(i-1)*p.N2Dp*p.N2Dp+4);
    end

    p.Kx = zeros(p.N2Dp,p.N2Dp*p.N_chi2D);
    for i = 1:p.N_chi2D
        p.Kx(1,(i-1)*p.N2Dp+1:(i-1)*p.N2Dp+2) = (i-1)*p.N2Dp*p.N2Dp+1:(i-1)*p.N2Dp*p.N2Dp+2;
        p.Kx(2,(i-1)*p.N2Dp+1:(i-1)*p.N2Dp+2) = fliplr((i-1)*p.N2Dp*p.N2Dp+3:(i-1)*p.N2Dp*p.N2Dp+4);
    end

    p.Kopplung_x = p.Kx(:,2:end-1);
    p.Kopplung_y = p.Ky(2:end-1,:);

    p.xi_K2D = linspace(-p.L_xi/2, p.L_xi/2, p.N_xi2D+1);
    for a = 1 : p.N_xi2D
        p.xi_Koord2D(1, p.N2Dp*(a-1)+1:p.N2Dp*a) = linspace(p.xi_K2D(a), p.xi_K2D(a+1), p.N2Dp);
    end
    p.dxi2D = diff(p.xi_K2D);
    p.dxi2D = p.dxi2D(1);
    p.chi_K2D = linspace(0, p.L_chi, p.N_chi2D+1);
    for a = 1 : p.N_chi2D
        p.chi_Koord2D(1, p.N2Dp*(a-1)+1:p.N2Dp*a) = linspace(p.chi_K2D(a), p.chi_K2D(a+1), p.N2Dp);
    end
    p.dchi2D = diff(p.chi_K2D);
    p.dchi2D = p.dchi2D(1);
    p.K_chi = meshgrid(p.chi_Koord2D,p.xi_Koord2D);
    p.K_xi = meshgrid(p.xi_Koord2D, p.chi_Koord2D)';
    [p.D_xi,p.phi, p.D] = eigendecomp2D(p);
    p.fermi2D = fermi2D(p.N_xi2D*p.N2Dp*p.N2Dp, p.phi);
end
end
