function [EfL, EfR] = solve_fermi(Esub, mat)



N_s = mat.N_s;
N_s_slice = N_s*1E+6*mat.W_c*1E-9;
N_d = mat.N_d; 
N_d_slice = N_d*1E+6*mat.W_c*1E-9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluation of the fermi-level at the left contact
% 1) Brute-Force for an rough initial Guess
EsubL = reshape(Esub(:,1,:), mat.n_of_valleys, mat.n_of_modes);

EfL = min(EsubL)-2 : 0.1 : min(EsubL)+2;
dC = zeros(1, length(EfL));


for IF = 1 : length(EfL)
    
    c0(IF)    = eval_slice_density(EsubL, EfL(IF), mat);
    dC(IF) = abs(log(c0(IF)-N_s_slice));
    
end

[~, idx] = min(dC);
EfL = EfL(idx);

% 2) Newton-Rhapson method for increased precision

newt_cond = 1E-4;
iter_max  = 250;
iter_err  = 1;
iter_cnt  = 0;

while iter_err > newt_cond && iter_cnt < iter_max 
    
    dc_dEf = eval_slice_density_derivative(EsubL, EfL, mat);
    c = eval_slice_density(EsubL, EfL, mat);
    EfL = EfL - (c-N_s_slice)/dc_dEf;
    
    iter_err = abs((c-N_s_slice)/N_s_slice);
    
    iter_cnt = iter_cnt+1;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluation of the fermi-level at the right contact
% 1) Brute-Force for an rough initial Guess
EsubR = reshape(Esub(:,mat.Nx,:), mat.n_of_valleys, mat.n_of_modes);

EfR = min(EsubR)-2 : 0.1 : min(EsubR)+2;
dC = zeros(1, length(EfL));


for IF = 1 : length(EfR)
    
    c0(IF) = eval_slice_density(EsubR, EfR(IF), mat);
    dC(IF) = abs(log(c0(IF)-N_d_slice));
    
end

[~, idx] = min(dC);
EfR = EfR(idx);

% 2) Newton-Rhapson method for increased precision

newt_cond = 1E-6;
iter_max  = 250;
iter_err  = 1;
iter_cnt  = 0;

while iter_err > newt_cond && iter_cnt < iter_max 
    
    dc_dEf = eval_slice_density_derivative(EsubR, EfR, mat);
    c = eval_slice_density(EsubR, EfR, mat);
    EfR = EfR - (c-N_d_slice)/dc_dEf;
    
    iter_err = abs((c-N_d_slice)/N_d_slice);
    
    iter_cnt = iter_cnt+1;
end
    
    



