function results = solve_transient_self_consistent_3D(mat, solver)

switch solver
    
    case 'SIGMA'
        
        fprintf('#####################################################\n');
        fprintf('#### TRANSIENT SIMULATION IS RUNNING ################\n');
        fprintf('Generation of the Initial Condition \n');
        
        Vd_ref = mat.Vd;
        Vg_ref = mat.Vg;
        
        mat.Vd = Vd_ref(1);
        mat.Vg = Vg_ref(1);
        
        fprintf('Initializing the biases: \t %.2d V at Gate and \t %.2d V at Drain', mat.Vg, mat.Vd);
        
        % calculation of the initial condition for the transient solver
        
        INIT = solve_self_consistent(mat, solver);
        % save('INIT.mat', 'INIT');
        
        %INIT = importdata('INIT.mat');
        mat.Vd = Vd_ref;
        mat.Vg = Vg_ref;
        
        if mat.CPL<2
            disp('tbd')
            results = solve_transient_SIGMA_SC_3D(mat, INIT);
        else
            disp('tbd')
            results = solve_transient_SIGMA_FC_3D(mat, INIT);
        end


    case 'DG'
        Temp = mat.Vg;
        mat.Vg = mat.Vg_INIT;    % set pre-switching gate voltage

        INIT = solve_self_consistent_3D(mat, solver);

        save('INIT.mat', 'INIT');

        % INIT = importdata('INIT.mat');
        
        mat.Vg = Temp;

        results = solve_transport_DG_3D_transient_SC_new(mat, INIT);

        
end
end