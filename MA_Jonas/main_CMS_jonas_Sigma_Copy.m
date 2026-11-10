% function main_CMS()
clear all; 
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
device = 'TAPSQDGFET'; %DGFET-type-I
solver = 'SIGMA'; %   <- hier zwischen NEGF_CMS_ALT, NEGF, und SIGMA wählen
simulation = 'self-consistent';% 'transient', 'self-consistent'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVICE PARAMETERS
dx      = 0.25; %0,25
dy      = 0.25; %0,25
w_ch    = 4;
offsetL = 1.0;
offsetU = 0.0;
taper   =   0; 
delete(gcp('nocreate'));
% parpool(8);
% parpool(6);
warning('off')
mat = buildDeviceCMS(device, dx,dy,w_ch,0.0,offsetL, offsetU,taper); %build
warning('on')

% SIMULATION PARAMETERS
mat.CPL             = 1;   % 0/1/2   für sigma auf 1
mat.n_of_modes      = 2;   % 1,2,3,...
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
% applied external biases 
switch simulation
    case 'self-consistent' 
        % Stationary Solver
        
        mat.Vs = +0.0;
%         mat.Vd = +0.1;
        mat.Vd = 0.1;
%         mat.Vd = fliplr(mat.Vd)
%         mat.Vd = [0.0 0.1 0.2 0.3]
        mat.Vg = +1.4;
        % mat.deltav = 0;0.325/2;
    case 'transient'
        % Transient Solver
        
        mat.t_max = 15E-12;
        mat.dt    = 2.5E-15;
        mat.t     = 0 : mat.dt : mat.t_max;
        %mat.Nt    = length(mat.t);
        mat.Vs    = +0.0;

        mat.Nt    = length(mat.t);
        mat.Vg    = 1.5311+0.336/2;
        mat.deltav = 0.0;
        mat.Vd    = 0.24+5e9*mat.t;
        
end
mat.transv_np = 0;
mat.FMC = 1e6;
mat.W_EVEN = 0;
mat.W_ODD = 0;
% poisson options

mat.poisson.opt.iter_err =  1E-3;
mat.poisson.opt.iter_max = 1;
mat.poisson.opt.alpha = 1;
mat.poisson.opt.solve = 'newton-rhapson';% 'newton-rhapson'; 'direct';

mat.wigner.params.Nk       = 200;
mat.wigner.params.k_max    = 4.5E+9;
mat.wigner.params.Ly       = pi/mat.wigner.params.k_max*(mat.wigner.params.Nk-1);
%pi/mat.wigner.params.k_max*(mat.wigner.params.Nk-1)/1E-9
mat.wigner.params.cap_ampl = 4.0;
mat.wigner.params.cap_rel  = 0.3;
mat.wigner.params.cap_n    = 4;

mat.sigma.params.Basis     = 'EXP'; %Ny = N+1 for continuity
mat.sigma.params.N         = 80; %100 working %Ny = N+1 for continuity
mat.sigma.params.Ny        = 601;%400;
mat.sigma.params.Ly        = 160e-9;
mat.sigma.params.cap_ampl  = 3.0; % 2 working
mat.sigma.params.cap_rel   = 0.2;
mat.sigma.params.cap_n     = 3;


mat.dE=2*2e-3;

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


if (strcmpi(solver,'SIGMA')==true)
    k = -2*pi/mat.sigma.params.Ly*((1:mat.sigma.params.N)-1/2*(mat.sigma.params.N+1));
end



if (strcmpi(solver,'NEGF_CMS_ALT')==true) || (strcmpi(solver,'NEGF')==true)
    RHO_Up = zeros((1/2)*size(results.RHO));
    RHO_Do = RHO_Up;
    RHO_Up_Do = RHO_Up;
    RHO_Do_Up = RHO_Up;
    
    
    for IN = 1:size(RHO_Do,1)
        for IM = 1:size(RHO_Do,1)
            RHO_Up(IN,IM)    = results.RHO(1+2*(IN-1), 1+2*(IM-1));
            RHO_Do(IN,IM)    = results.RHO(2*IN, 2*IM);
            RHO_Up_Do(IN,IM) = results.RHO(2*IN, 1+2*(IM-1));
            RHO_Do_Up(IN,IM) = results.RHO(1+2*(IN-1), 2*IN);
        end
    end
    
    
    
    G = sum(results.G, 3);
    
    
    G_Up = zeros((1/2)*size(G));
    G_Do = G_Up;
    G_Up_Do = G_Up;
    G_Do_Up = G_Up;
    
    
    
    for IN = 1:size(G_Up,1)
        for IM = 1:size(G_Up,1)
            G_Up(IN,IM) = G(1+2*(IN-1), 1+2*(IM-1));
            G_Do(IN,IM) = G(2*IN, 2*IM);
        end
    end
    
    
    %Diag = diag(flipud(abs(results.RHO(1:2:end,1:2:end))),0);
    %Diag1 = diag(flipud(abs(RHO_Up)),30);
    
    %[RHO_Up_RS, r_up, r_strich_up] = Schwerpunkttrafo(size(RHO_Do,1), RHO_Up);    %Anders lösen -> xis und chis einzeln berechnen und nur einen Vektor rausnehmen
    %[RHO_Do_RS, r_do, r_strich_do] = Schwerpunkttrafo(size(RHO_Do,1), RHO_Do);
    
%% build sparse structure to get transformed coordinates  [x', x, xii, chii, rho_up, rho_down], [rho_up, rho_down, rho_up_down, rho_down_up]

    dE = mat.dE;
    counter = size(RHO_Do,1);

    Frame = zeros(counter*counter,4);
    Frame_rho_values = zeros(counter*counter,4);
    row_vec=zeros(counter,1);


    for i=1:counter
        row_vec(i,1)=i;
    end


    for i=1:counter
        Frame(((i-1)*counter+1):i*counter,1) = i*ones(counter,1);
        Frame(((i-1)*counter+1):i*counter,2) = row_vec;
    end


   for i=1:(counter*counter)
       Frame_rho_values(i,1) = RHO_Up(Frame(i,1),Frame(i,2));
       Frame_rho_values(i,2) = RHO_Do(Frame(i,1),Frame(i,2));
       Frame_rho_values(i,3) = RHO_Up_Do(Frame(i,1),Frame(i,2));
       Frame_rho_values(i,4) = RHO_Do_Up(Frame(i,1),Frame(i,2));
       Frame(i,3) = Frame(i,2) - Frame(i,1);
       Frame(i,4) = (Frame(i,2) + Frame(i,1))/2;
   end

    
   searcher=find(Frame(:,4)==19);
   rho_slice_up    = Frame_rho_values(searcher,1);
   rho_slice_do    = Frame_rho_values(searcher,2);
   rho_slice_up_do = Frame_rho_values(searcher,3);
   rho_slice_do_up = Frame_rho_values(searcher,4);


   fill_up = counter - size(rho_slice_up,1);

   full_slice_up    = [zeros(fill_up/2,1);rho_slice_up;zeros(fill_up/2,1)];
   full_slice_do    = [zeros(fill_up/2,1);rho_slice_do;zeros(fill_up/2,1)];
   full_slice_up_do = [zeros(fill_up/2,1);rho_slice_up_do;zeros(fill_up/2,1)];
   full_slice_do_up = [zeros(fill_up/2,1);rho_slice_do_up;zeros(fill_up/2,1)];
   
  %% Q-Matrix for Fourier Transformation
   
   Ny_alt = mat.sigma.params.Ny;
   Ny     = size(RHO_Do,1);
   Nk     = 80;
   Ly     = mat.sigma.params.Ly;
   y      = linspace(-Ly/2, +Ly/2, Ny);
    
   Q = zeros(Ny, Nk);
            
   k     = -2*pi/Ly*((1:Nk)-1/2*(Nk+1));
   dk    = abs(k(2)-k(1));
    
   for IN = 1 : Nk
       Q(:,IN) = 1./sqrt(Ly)*exp(+1i*k(IN)*y');
   end
    
   wf_slice_up = Q'*full_slice_up;
   wf_slice_do = Q'*full_slice_do;

  % %% Q Alternative
  % 
  %  Ny_alt = mat.sigma.params.Ny;
  %  %Ny     = size(RHO_Do,1);
  %  Nk     = 80;
  %  Ly     = mat.sigma.params.Ly;
  %  y_alt      = linspace(-Ly/2, +Ly/2, Ny_alt);
  % 
  %  Q_alt = zeros(Ny_alt, Nk);
  % 
  %  k_alt     = -2*pi/Ly*((1:Nk)-1/2*(Nk+1));
  %  dk_alt    = abs(k_alt(2)-k_alt(1));
  % 
  %  for IN = 1 : Nk
  %      Q_alt(:,IN) = 1./sqrt(Ly)*exp(+1i*k_alt(IN)*y_alt');
  %  end
  % 
  % 
  %  diff = Ny_alt - counter;
  %  full_slice_up_alt = [zeros((fill_up+diff)/2,1);rho_slice_up;zeros((fill_up+diff)/2,1)];
  % 



   %f_Up = Q'*RHO_Up_RS;
   %f_Do = Q'*RHO_Do_RS;
    
   %figure(1), mesh(abs(f_Up));
   %figure(2), mesh(abs(f_Do));

end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if mat.saveres
%     save(strcat('/work/smmapech/msascatter/res/STAT_dev7nm_SQ_',int2str(mat.n_of_modes),'nm_',num2str(mat.enableSD),'_SD',num2str(mat.K_order),'k_order',num2str(mat.sigma.params.N.mat),'N.mat'),'results')
% end
