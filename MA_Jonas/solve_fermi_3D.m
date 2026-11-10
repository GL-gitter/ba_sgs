function [EfL, EfR] = solve_fermi_3D(Esub, mat)

%% Ändern

N_s = mat.N_s;
N_s_slice = N_s*1E+6*mat.W_cy*1E-9*mat.W_cz*1E-9;
N_d = mat.N_d; 
N_d_slice = N_d*1E+6*mat.W_cy*1E-9*mat.W_cz*1E-9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluation of the fermi-level at the left contact
% 1) Brute-Force for an rough initial Guess
EsubL = reshape(Esub(:,1,:), mat.n_of_valleys, mat.n_of_modes);

% EfL = min(EsubL)-2 : 0.1 : min(EsubL)+2;
EfL = linspace(0, min(EsubL)+2, 201);
dC = zeros(1, length(EfL));


for IF = 1 : length(EfL)
    
    c0(IF)    = eval_slice_density_3D(EsubL, EfL(IF), mat);
    dC(IF) = abs(log(c0(IF)-N_s_slice));
    
end

[~, idx] = min(dC);
EfL = EfL(idx);

% 2) continue brute force
rbf = 0.1;
for IBF =1:5    
    EfL = linspace(EfL-rbf,EfL+rbf, 201);
    dC = zeros(1, length(EfL));
    for IF = 1 : length(EfL)
        c0(IF)    = eval_slice_density_3D(EsubL, EfL(IF), mat);
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

EfR = linspace(0, min(EsubR)+2, 201);
dC = zeros(1, length(EfR));


for IF = 1 : length(EfR)
    
    c0(IF)    = eval_slice_density_3D(EsubR, EfR(IF), mat);
    dC(IF) = abs(log(c0(IF)-N_s_slice));
    
end

[~, idx] = min(dC);
EfR = EfR(idx);

% 2) continue brute force
rbf = 0.1;
for IBF =1:5    
    EfR = linspace(EfR-rbf,EfR+rbf, 201);
    dC = zeros(1, length(EfR));
    for IF = 1 : length(EfR)
        c0(IF)    = eval_slice_density_3D(EsubR, EfR(IF), mat);
        dC(IF) = abs(log(c0(IF)-N_s_slice));
    end
    
    [~, idx] = min(dC);
    EfR = EfR(idx);
    rbf=rbf*0.1;
end    
    



