function [A,rhs] = get_SysM(mat, p, EM)

[A,rhs] = get_Diff(p);


if (strcmpi(mat.dg.params.Xi_solver,'FV')==true)
    G_glob = get_Drift_accel_trial_permutedShell(mat, p,0, EM);       %change to get_Drift_accel_trial
elseif (strcmpi(mat.dg.params.Xi_solver,'DG')==true)
    if (strcmpi(mat.dg.params.Q,'EIG')==true)
      G_glob = conj(get_Drift_accel_trial_permutedShell(mat, p,0, EM));
    elseif (strcmpi(mat.dg.params.Q,'EXP')==true)
      G_glob = get_Drift_accel_trial_permutedShell(mat, p,0, EM); 
    end
end

A = A+G_glob;


%% Plot Flux

plotFlux = false;

if plotFlux == true
    plotColors = cell(1,4);
    
    plotColors{1} = 'green';
    plotColors{2} = 'blue';
    plotColors{3} = 'red';
    plotColors{4} = 'yellow';
    
    
    close 
    

    %%%% Plot Upwind left to right
    fluxDiagUp = zeros(p.N_K_chi, p.Nk);
    
    blockOffset = 20;
    
    for k=1:p.N_K_chi
        for i=1:p.Nk
            fluxDiagUp(k,i)=A((blockOffset-1)*(p.N_K_chi*p.Nk)+(k-1)*p.Nk+i, blockOffset*(p.N_K_chi*p.Nk)+i);
        end
    end
        
    
    fluxFigUp = figure('Name','FluxUp');
    
    figure(fluxFigUp)
    for k=1:p.N_K_chi
        plot(squeeze(fluxDiagUp(k,:)),'LineWidth',2,'Color',plotColors{k});
        hold on
    end
    
    
    
    %%%% Plot Downwind right to left
    
    fluxDiagDown = zeros(p.N_K_chi, p.Nk);
    
    blockOffset = 20;
    
    for k=1:p.N_K_chi
        for i=1:p.Nk
            fluxDiagDown(k,i)=A(blockOffset*(p.N_K_chi*p.Nk)+i+(k-1)*p.Nk , (blockOffset-1)*(p.N_K_chi*p.Nk)+(p.N_K_chi-1)*(p.Nk)+i);
        end
    end
    
    fluxFigDown = figure('Name','FluxDown');
    
    figure(fluxFigDown)
    for k=1:p.N_K_chi
        plot(squeeze(fluxDiagDown(k,:)),'LineWidth',2,'Color',plotColors{k});
        hold on
    end
end

%%



end


