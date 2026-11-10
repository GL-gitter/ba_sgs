function [EfL, EfR, mexv] = solve_fermi_3D_MP(Esub, Vsub, mat)



N_s = mat.N_s;
N_s_slice = N_s*1E+6*mat.W_cy*1E-9*mat.W_cz*1E-9;
N_d = mat.N_d; 
N_d_slice = N_d*1E+6*mat.W_cy*1E-9*mat.W_cz*1E-9;


includemefffermi = 1;
if includemefffermi == 1

    h    = 6.626E-34;
    m0   = 9.109E-31;
    hb   = h/2/pi;
    q    = 1.602E-19;
    dx  = mat.dx*1e-9;
    Ny=mat.Ny;
    Nx=mat.Nx;
    Nz=mat.Nz;

    t0 = hb^2/2/m0/dx^2/q;

    H     = eval_full_Hamiltonian_MV_alt_3D(mat,mat.V*0);

    VMS = sparse(Nx*Ny*Nz,Nx*mat.n_of_modes);

    for IX=0:Nx-1                                               %this is kinda stupid but reshaping in matlab sucks
        vtemp=squeeze(Vsub(1,IX+1,:,:,:));
        vtemp2=[];
        for IM=1:mat.n_of_modes
            vtemp3=squeeze(vtemp(IM,:,:))';
            vtemp2=[vtemp2,vtemp3(:)];
        end
        VMS(IX*Ny*Nz+1:(IX+1)*Ny*Nz,IX*mat.n_of_modes+1:(IX+1)*mat.n_of_modes) = vtemp2;
    end


    HMS_L = squeeze(VMS(1:Ny*Nz,1:mat.n_of_modes))' * H(1:Ny*Nz,Ny*Nz+1:2*Ny*Nz) *squeeze(VMS(1:Ny*Nz,1:mat.n_of_modes));
    % HMS_R = squeeze(VMS(1,end,:,:))' *H(end-Ny+1:end,end-2*Ny+1:end-Ny)* squeeze(VMS(1,end,:,:));
    mexv=-t0./diag(HMS_L);
    % tR=-t0./diag(HMS_R);

else
    mexv = mat.me_x_ch*ones(mat.n_of_modes,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluation of the fermi-level at the left contact
% 1) Brute-Force for an rough initial Guess
EsubL = reshape(Esub(:,1,:), mat.n_of_valleys, mat.n_of_modes);

% EfL = min(EsubL)-2 : 0.1 : min(EsubL)+2;
EfL = linspace(0, min(EsubL)+2, mat.Nx);
dC = zeros(1, length(EfL));


for IF = 1 : length(EfL)
    
    c0(IF)    = eval_slice_density_3D_MP(EsubL,mexv, EfL(IF), mat);
    dC(IF) = abs(log(c0(IF)-N_s_slice));
    
end

[~, idx] = min(dC);
EfL = EfL(idx);

% 2) continue brute force
rbf = 0.1;
for IBF =1:5    
    EfL = linspace(EfL-rbf,EfL+rbf, mat.Nx);
    dC = zeros(1, length(EfL));
    for IF = 1 : length(EfL)
        c0(IF)    = eval_slice_density_3D_MP(EsubL,mexv, EfL(IF), mat);
        dC(IF) = abs(log(c0(IF)-N_s_slice));
    end
    
    [~, idx] = min(dC);
    EfL = EfL(idx);
    rbf=rbf*0.1;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluation of the fermi-level at the right contact
% 1) Brute-Force for an rough initial Guess
EsubR = reshape(Esub(:,mat.Nx,:), mat.n_of_valleys, mat.n_of_modes);

EfR = linspace(0, min(EsubR)+2, mat.Nx);
dC = zeros(1, length(EfR));


for IF = 1 : length(EfR)
    
    c0(IF)    = eval_slice_density_3D_MP(EsubR, mexv, EfR(IF), mat);
    dC(IF) = abs(log(c0(IF)-N_s_slice));
    
end

[~, idx] = min(dC);
EfR = EfR(idx);

% 2) continue brute force
rbf = 0.1;
for IBF =1:5    
    EfR = linspace(EfR-rbf,EfR+rbf, mat.Nx);
    dC = zeros(1, length(EfR));
    for IF = 1 : length(EfR)
        c0(IF)    = eval_slice_density_3D_MP(EsubR,mexv, EfR(IF), mat);
        dC(IF) = abs(log(c0(IF)-N_s_slice));
    end
    
    [~, idx] = min(dC);
    EfR = EfR(idx);
    rbf=rbf*0.1;
end    
    



