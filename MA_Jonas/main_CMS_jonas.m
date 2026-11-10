% function main_CMS()
addpath ./DG_code/
addpath ./DG_code/core/

clear workspace; 
% close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
device = 'DGFET-type-I'; %DGFET-type-I (uncoupled MS), DGFET (self/full coupled MS), GAAFET (3D MSA)
solver = 'DG'; %   <- hier zwischen NEGF_CMS_ALT, NEGF, SIGMA, WIGNER und DG wählen
simulation = 'self-consistent';% 'transient', 'self-consistent'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%mat.n_of_modes      = 3;   % 1,2,3,...
% DEVICE PARAMETERS
dx      = 0.25; %0,25 %Mind the following condition for DG: (Nx-1)/2 has to be even
dy      = 0.25; %0,25
dz      = 0.25;

%w_ch    = 5;
w_chy    = 3;
w_chz    = 3;
w_oxy    = 1;
w_oxz = w_oxy;


offsetL = 1.0;
offsetU = 0.0;
taper   =   0; 

delete(gcp('nocreate'));
% parpool(8);
% parpool(6);
warning('off')

%% Set up the device to be simulated

if (strcmpi(device,'GAAFET')==true)
    mat = buildDeviceCMS_3D(device, dx, dy, dz, w_chy, w_oxy, w_chz, w_oxz);

elseif (strcmpi(device,'DGFET')==true)
    mat = buildDeviceCMS(device, dx,dy,w_ch,0.0,offsetL, offsetU,taper);
elseif (strcmpi(device,'DGFET-type-I')==true)
    mat = buildDevice(device);
end

warning('on')

%%

%% SIMULATION PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat.MP              = 0;   % True if simulating for Mathias Paper
if mat.MP == 1
    disp('Running code for Mathias Paper!')
    disp('Permittivity increased towards contacts')
    lenp = 0.1; % in percent
    fakp = 100;
    lenp = ceil(lenp*mat.Nx);
    fil = ((fakp-1)*0.5*(1+cos((0:(lenp-1))/(lenp-1)*pi))+1);
    % plot(fil)
    mat.Eps(1:lenp,:,:) = mat.Eps(1:lenp,:,:) .* fil';
    mat.Eps(end-lenp+1:end,:,:) = mat.Eps(end-lenp+1:end,:,:) .* fliplr(fil)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat.CPL             = 0;   % 0/1/2   für sigma auf 1
mat.n_of_modes      = 1;   % 1,2,3,...
mat.doubleNm        = 0;   % 0/1
mat.toggleINIT      = 0;   % 0/1
mat.enableSD        = 0;   % 0/1
mat.VNBC            = 1;   % 0/1
mat.K_order         = 1;   % 1/2
mat.plot            = 0;   % 0/1
mat.saveres         = 0;   % 0/1
mat.time_max        = 0;   % t in h or 0 for inf
mat.Hamiltonian2D   = 0;   % 0/0.5/1
mat.reorder         = 0;
mat.smoothEM        = 0;
mat.smoothVM        = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% applied external biases 
switch simulation
    case 'self-consistent' 
        % Stationary Solver
        
        mat.Vs = +0.0;
%         mat.Vd = +0.1;
        mat.Vd = 1.0; %+1.0;  %[(0.0:0.05:+0.5) (0.55:0.15:1.3)]; %(0.0:0.05:+1.2);
%         mat.Vd = fliplr(mat.Vd)
%         mat.Vd = [0.0 0.1 0.2 0.3]
        mat.Vg = +1.4; %(+0.8:0.05:+1.6); %(0.8:0.05:+1.6); %[+1.2 +1.3 +1.4];
        % mat.deltav = 0;0.325/2;
    case 'transient'
        % Transient Solver
        
        mat.t_max = 3.0E-13;
        mat.dt    = 2E-17;
        mat.t     = 0 : mat.dt : mat.t_max;
        %mat.Nt    = length(mat.t);
        mat.Vs    = +0.0;

        mat.Nt      = length(mat.t);
        mat.Vg_INIT = 1.2;
        mat.Vg      = 1.6*ones(1,mat.Nt);
        Vg_temp     = linspace(1.2,1.6,50);
        mat.Vg(1:length(Vg_temp)) = Vg_temp;
        mat.deltav  = 0.0;
        mat.Vd      = 1.0;
        
end
mat.transv_np = 0;
mat.FMC = 1e6;
mat.W_EVEN = 0;
mat.W_ODD = 0;
% poisson options

mat.poisson.opt.iter_err =  1.0E-3;
mat.poisson.opt.iter_max = 25;
mat.poisson.opt.alpha = 1;
mat.poisson.opt.solve = 'newton-rhapson';% 'newton-rhapson'; 'direct';

%% Ab hier alles so lassen

%% Solver specific parameters

mat.wigner.params.Nk       = 100;   %200
mat.wigner.params.k_max    = 4.5E+9;
mat.wigner.params.Ly       = pi/mat.wigner.params.k_max*(mat.wigner.params.Nk-1);
%pi/mat.wigner.params.k_max*(mat.wigner.params.Nk-1)/1E-9
mat.wigner.params.cap_ampl = 4.0;
mat.wigner.params.cap_rel  = 0.3;
mat.wigner.params.cap_n    = 4;

mat.sigma.params.Basis     = 'EXP'; %Ny = N+1 for continuity
mat.sigma.params.N         = 80; %100 working %Ny = N+1 for continuity    % 80 for good current density results
mat.sigma.params.Ny        = 401;%400;                                    % 401 for good current density results
mat.sigma.params.Ly        = 160e-9;
mat.sigma.params.cap_ampl  = 3.0; % 2 working
mat.sigma.params.cap_rel   = 0.2;
mat.sigma.params.cap_n     = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DG parameters
if (strcmpi(simulation,'transient')==true)
    mat.dg.params.refresh_init_trans = input('Refresh of initial Data for transient Simulation? ');
else
    mat.dg.params.refresh_init_trans = false;
end


mat.dg.params.permute_shell = true;
mat.dg.params.rho = true;                                   % false: phase/eigen basis, true: relative-coordinate rho basis (real space)
mat.dg.params.rho_flux = 'rusanov';         % rho solver: 'central', legacy 'matrix-upwind', 'rusanov' or 'rusanov-boundary-upwind'
mat.dg.params.rho_diff_scale = 1;
mat.dg.params.rho_drift_scale = 1e9;

mat.dg.params.N_K_chi        = 3;    % Don't change
mat.dg.params.sl             = 0;

mat.dg.params.Xi_solver      = 'FV';    % FV ord DG

if (strcmpi(simulation,'transient')==true)
    mat.dg.params.st_init    = true;
else
    mat.dg.params.st_init    = false;
end

if mat.dg.params.rho
    mat.dg.params.N_xi           = 300;   %160 for real space based approach gives best results
else
    mat.dg.params.N_xi           = 400;   %400;   %Choose even number for DG!
end

if (strcmpi(mat.dg.params.Xi_solver,'DG')==true)
    mat.dg.params.N_xi       = 80;                %133 for 'EXP'
end
mat.dg.params.Q              = 'EXP'; % 'EXP or 'EIG'
mat.dg.params.N_K_xi         = 3;     % default is 3
mat.dg.params.Ly             = 160e-9;

mat.dg.params.N              = 80;   % Number of elements in k-direction

mat.dg.params.transient_mode = 'self-consistent';     % 'flatband' or 'self-consistent', flatband just for development to test dimensions and stuff like that
if (strcmpi(simulation,'transient')==true)
    mat.dg.params.MPC            = 2.0e-10;    % (2.0e-12) for transient simulation
elseif (strcmpi(simulation,'self-consistent')==true)
    mat.dg.params.MPC            = 1;
end

if ((strcmpi(mat.dg.params.Q,'EXP')) && (strcmpi(mat.dg.params.Xi_solver,'DG'))==true)
    mat.dg.params.Nk_K       = 4;      % Nk_K is number of nodes per element for xi-direction when using exp base for fourier transform
else
    mat.dg.params.Nk_K       = 1;
end

%mat.dg.params.cap_ampl  = 3.0; % 2 working
%mat.dg.params.cap_rel   = 0.2;
%mat.dg.params.cap_n     = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
mat.dE=2*2e-3;        %% Changed for Mathias Paper; default is 2*2e-3

mat.Tau = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization of the log-file

% fileID = fopen(['log_file_', device, '_', solver, '_', simulation, '.txt'], 'w');

%initialize_log_file(mat, fileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try mat.Hamiltonian2D;
catch
    mat.Hamiltonian2D=0;
end

try mat.CPL;
catch
    mat.CPL=0;
end

try mat.CPL;
catch
    mat.CPL=0;
end

try mat.doubleNm;
catch
    mat.doubleNm=0;
end

try mat.debug;
catch
    mat.debug=0;
end



if mod(mat.sigma.params.Ny,2)==1 && strcmpi(mat.sigma.params.Basis,'EIG')==true
    disp('Use even amounts of points in xi with eigen basis!')
    mat.sigma.params.Ny = mat.sigma.params.Ny+1;
end

results = solve(mat, solver, simulation);
results.mat = mat;

save('results.mat','results');


%figure(444),mesh(results.chi,results.xi,real(squeeze(results.rho(1,1,:,:)),FaceColor='interp')



