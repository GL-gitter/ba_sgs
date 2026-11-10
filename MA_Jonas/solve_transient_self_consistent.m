function results = solve_transient_self_consistent(mat, solver)

switch solver
    
    case 'NEGF'
        
        fprintf('NEGF is no valid solver for the transient simulation\n');
        error('NEGF is no valid solver for the transient case\n');
        
    case 'WIGNER'
        
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
        
        % INIT = importdata('INIT.mat');
        mat.Vd = Vd_ref;
        mat.Vg = Vg_ref;
        
        results = solve_transient_WIGNER_accel(mat, INIT);

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
        
        % INIT = importdata('INIT.mat');
        mat.Vd = Vd_ref;
        mat.Vg = Vg_ref;
        if mat.CPL<2
            results = solve_transient_SIGMA_SC(mat, INIT);
        else
            if mat.doubleNm==0
                results = solve_transient_SIGMA_FC(mat, INIT);
            else
                results = solve_transient_SIGMA_FC_DNM(mat, INIT);
            end
        end
        
    case 'OMICRON'
        
        fprintf('#####################################################\n');
        fprintf('#### TRANSIENT SIMULATION IS RUNNING ################\n');
        fprintf('Generation of the Initial Condition \n');
        
        Vd_ref = mat.Vd;
        Vg_ref = mat.Vg;
        
        mat.Vd = Vd_ref(1);
        mat.Vg = Vg_ref(1);
        
        fprintf('Initializing the biases: \t %.2d V at Gate and \t %.2d V at Drain', mat.Vg, mat.Vd);
        
        % calculation of the initial condition for the transient solver
        
        INIT = solve_self_consistent_OMICRON(mat, solver);
        % save('INIT.mat', 'INIT');
        
        % INIT = importdata('INIT.mat');
        mat.Vd = Vd_ref;
        mat.Vg = Vg_ref;

        results = solve_transient_OMICRON_accel(mat, INIT); 


    case 'SIGMA_SCATTER'
        
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
       %save('INIT.mat', 'INIT');
        
        %INIT = importdata('INIT.mat');
        mat.Vd = Vd_ref;
        mat.Vg = Vg_ref;
        
        results = solve_transient_S_SIGMA_accel(mat, INIT);

    case 'DG'
        if (strcmpi(mat.dg.params.transient_mode,'flatband')==true)
            %% Flatband: get stationary initial data first (system matrix, boundary condition) befor switching
            Temp = mat;
            Vg = mat.Vg(1);    % set pre-switching gate voltage

            INIT = get_init_DG_2D_trans_flatband(mat, Vg);
            
            mat = Temp;
            results = solve_transport_DG_2D_transient_flatband(mat, INIT);

        elseif (strcmpi(mat.dg.params.transient_mode,'self-consistent')==true)
            Temp = mat.Vg;

            if mat.dg.params.refresh_init_trans
                mat.Vg = mat.Vg_INIT;    % set pre-switching gate voltage
    
                INIT = solve_self_consistent(mat, solver);
    
                save('INIT.mat', 'INIT');
            else
                INIT = importdata('INIT.mat');
            end
            
            mat.Vg = Temp;

            if isfield(mat.dg.params, 'rho') && mat.dg.params.rho
                results = solve_transport_DG_2D_transient_SC_new_rho(mat, INIT);
            else
                results = solve_transport_DG_2D_transient_SC_new(mat, INIT);
            end
            
        end
        
end
