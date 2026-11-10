function results = solve(mat, solver, simulation)

switch simulation
    
    case 'self-consistent'
        if strcmp(solver, 'OMEGA') || strcmp(solver, 'OMICRON') 
            results = solve_self_consistent_OMICRON(mat, solver);   
        else
            if strcmp(mat.type, 'GAAFET')
                results = solve_self_consistent_3D(mat, solver);  
            else
                if mat.time_max==0
                    results = solve_self_consistent(mat, solver);  
                else
                    results = solve_self_consistent_ls(mat, solver);  
                end
            end
        end
        
    case 'flatband'
        results = solve_flatband(mat, solver);
        
    case 'transient'
        if strcmp(mat.type, 'GAAFET')
            results = solve_transient_self_consistent_3D(mat, solver);
        else
            results = solve_transient_self_consistent(mat, solver);
        end
        
    otherwise
        error('no valid solver determined\n');

end