function results = solve_self_consistent_3D(mat, solver)

%% If scattering is enabled -> call itself for equilibrium distribution
% if mat.Tau ~= 0
%     fprintf('Get Equilibrium distribution');
%     tautemp = mat.Tau;
%     vdtemp  = mat.Vd;
%     mat.Vd = 0;
%     mat.Tau = 0;
%     % vgtemp = mat.Vg;
% 
%     resultseq = solve_self_consistent(mat, solver);   
% 
%     mat.f_eq=zeros(mat.n_of_modes, mat.sigma.params.N, length(mat.x));
%     for IM =1:mat.n_of_modes
%         mat.f_eq(IM,:,:) = squeeze(resultseq.sigT(IM,:,:))./resultseq.nm(IM,:);
%     end
%     results.feq = mat.f_eq;
%     mat.Tau = tautemp;
%     mat.Vd  = vdtemp;
%     clear resultseq;
% end

%% Adaption for DG algorithm
N_K_chi = mat.dg.params.N_K_chi;
N_chi   = ((mat.Nx-1)/(N_K_chi-1));
mat.Np      = N_chi*N_K_chi;
%%

% solve for the subbands

[Esub, Vsub] = solve_subbands_3D(mat);

% solve for the fermi level

if mat.MP == 1
    [Ef0_L, Ef0_R, mexv] = solve_fermi_3D_MP(Esub, Vsub, mat);
    mat.mexv = mexv;
else
    [Ef0_L, Ef0_R] = solve_fermi_3D(Esub, mat);
end
% if mat.Tau ~= 0
%     mat.Vg=mat.deltav;
% else
%     mat.Vg=0;
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SELF CONSISTENT LOOP

[M, node_num, node_type] = eval_Poisson_3D(mat);


if (strcmpi(solver,'DG')==true)   %% necessary cause DG results are not structured in spatial discretization (two values at cell borders)
    n_total = zeros(mat.Np, mat.Ny, mat.Nz, length(mat.Vd), length(mat.Vg));
    j_total = zeros(mat.Np, length(mat.Vd), length(mat.Vg));
    V_total = zeros(mat.Nx, mat.Ny, mat.Nz, length(mat.Vd), length(mat.Vg));
else
    n_total = zeros(mat.Nx, mat.Ny, mat.Nz, length(mat.Vd),length(mat.Vg));
    j_total = zeros(mat.Nx, length(mat.Vd),length(mat.Vg));
    V_total = zeros(mat.Nx, mat.Ny, mat.Nz, length(mat.Vd), length(mat.Vg));
% jm = 0;
end


performance_switch = 0;

for IG = 1 : length(mat.Vg)
    
    
    Vg_bias = Ef0_L + mat.phi_m_g-mat.Xi-mat.Vg(IG); % -mat.Vg(IG); % Matched work functions  v + mat.phi_m_g-mat.Xi-mat.Vg(IG);
    
%     fprintf(fileID, '###########################################################\n');
%     fprintf(fileID, 'Applied Gate Voltage:\t %.2f V \n', mat.Vg(IG));
    % CPL=2;
    CPL=0;
    for ID = 1 : length(mat.Vd)
        
        if ID > 2
            
            Vbi = 2*V_total(:,:,:,ID-1)-V_total(:,:,:,ID-2);
            
%             mat.poisson.opt.iter_err = 1E-4;
%             mat.poisson.opt.iter_max = 15;
        else
            Vbi=zeros(mat.Nx,mat.Ny,mat.Nz);
            Vbi(mat.gate_idx,:,:)=Vg_bias;
            Vbi(mat.x>=mat.L_s+mat.L_c,:,:) = -mat.Vd(ID);   
        end
        
        iter_cnt = 0;
        iter_err = 1;
        
        while iter_err > mat.poisson.opt.iter_err && iter_cnt < mat.poisson.opt.iter_max
            
            % if length(mat.Vd)==1
            % 
            %     if iter_cnt==1 & mat.CPL>0
            %         CPL=1;
            %         disp('Self coupling turned on');
            %     end
            % 
            %     if iter_cnt==3 & mat.CPL==2
            %         CPL=2;
            %         disp('Full coupling turned on');
            %     end
            % 
            % else
            %     if ID < 3 & mat.CPL>0
            %         CPL=1;
            %         disp('Self coupling turned on');
            %     elseif ID==3 & mat.CPL>1
            %         CPL=2;
            %         disp('Full coupling turned on');
            %     end
            % end

            Vref = mat.V; mat.V = Vref+Vbi;
            [EM, VM] = solve_subbands_3D(mat);
            
            %K = getK_3D(mat.Nx,mat.dx*1e-9,mat.n_of_modes, 1,mat.me_x_ch, mat.me_x, VM, mat.K_order);
            % 
            % if mat.plot; plotK(K,44,mat.n_of_modes); end
            % figure(33),plot(squeeze(EM))
            
            % disp('K0 set to 0');
            % for IMO=1:mat.n_of_modes
            %     for IMI=1:mat.n_of_modes
            %         if IMO ~= IMI
            %             mat.K0(IMO,IMI,:)=0;
            % 
            %         end
            %     end
            % end
            
            mat.V = Vref;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % MODE SPACE MAIN ALGORITHM
            
            switch solver
                
                % case 'NEGF'
                % 
                %     NEGF = solve_transport_NEGF(mat, EM, VM, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID));
                % 
                %     n = NEGF.n;
                %     j = NEGF.j;
                % 
                % case 'NEGF_RS'
                % 
                %     if (length(mat.Vd)>3 && ID<3) || (length(mat.Vd)==1 && iter_cnt < 10)
                %         disp('MS Speed up')
                %         NEGF = solve_transport_NEGF(mat, EM, VM, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID));
                %         n = NEGF.n;
                %         j = NEGF.j;
                % 
                %     else
                %         tic;
                %         NEGF = solve_transport_NEGF_real_space(mat, Vref+Vbi, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID));
                %         toc
                %         n = NEGF.n;
                %         j = NEGF.j;
                %     end
                % 
                % case 'NEGF_RSMV'
                % 
                %     if (length(mat.Vd)>3 && ID<3) || (length(mat.Vd)==1 && iter_cnt < 10)
                %         disp('MS Speed up')
                %         NEGF = solve_transport_NEGF(mat, EM, VM, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID));
                %         n = NEGF.n;
                %         j = NEGF.j;
                %     else
                %         tic;
                %         NEGF = solve_transport_NEGF_real_space_MV(mat, Vref+Vbi, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID));
                %         toc
                %         n = NEGF.n;
                %         j = NEGF.j;
                %     end
                % 
                 case 'WIGNER'
                 
                     if performance_switch == 0
                 
                 
                         WIGNER = solve_transport_WIGNERconst_eff_optim_3D(mat, EM, VM, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID));
                         results.C_wo_drift = WIGNER.C_wo_drift;
                         performance_switch = 0;
                 
                 
                     else
                         WIGNER = solve_transport_WIGNERconst_eff_accel(results.C_wo_drift,mat, EM, VM, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID));
                 
                     end
                 
                     n = WIGNER.n;
                     j = WIGNER.j;
                     results.wf = WIGNER.wf;
                     results.Ef0 = Ef0_L;
                
                case 'SIGMA'
                        if CPL == 0
                            clear('K');
                            K.K0 = ones(mat.n_of_modes,mat.n_of_modes,mat.Nx);
                            K.K1 = zeros(mat.n_of_modes,mat.n_of_modes,mat.Nx);
                            K.K2 = zeros(mat.n_of_modes,mat.n_of_modes,mat.Nx);
                            SIGMA = solve_transport_SIGMA_SC_3D(mat, EM, VM, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID),K);
                            n = SIGMA.n;
                            j = SIGMA.j;
                        
                        elseif CPL == 1
                            % K.K1=K.K1*0;
                            % K.K2=K.K2*0;
                            tic
                            SIGMA = solve_transport_SIGMA_SC_3D(mat, EM, VM, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID),K);
                            toc
                            n = SIGMA.n;
                            j = SIGMA.j;
                        
                        elseif CPL == 2
                            if mat.doubleNm == 1
                                if mat.n_of_modes==6
                                    vecp1 = [1,3,5];
                                    vecp2 = vecp1+1;
                                elseif mat.n_of_modes==5
                                    vecp1 = [1,3,5];    
                                    vecp2 = [2,4];   
                                elseif mat.n_of_modes==4
                                    vecp1 = [1,3];
                                    vecp2 = vecp1+1;
                               elseif mat.n_of_modes==3
                                    vecp1 = [1,3];
                                    vecp2 = [2];
                                end

                                EM_ODD  = EM(1,:,vecp1);
                                EM_EVEN = EM(1,:,vecp2);
        
                                VM_ODD  = VM(1,:,vecp1,:,:);
                                VM_EVEN = VM(1,:,vecp2,:,:);

                                
                                K_ODD.K0 = K.K0(vecp1,vecp1,:);
                                K_ODD.K1 = K.K1(vecp1,vecp1,:);
                                K_ODD.K2 = K.K2(vecp1,vecp1 ,:);

                                K_EVEN.K0 = K.K0(vecp2,vecp2,:);
                                K_EVEN.K1 = K.K1(vecp2,vecp2,:);
                                K_EVEN.K2 = K.K2(vecp2,vecp2,:);

                                mat.n_of_modes=ceil(length(vecp2));

                                if mat.n_of_modes>1
                                    SIGMA_EVEN = solve_transport_SIGMA_FC_3D(mat, EM_EVEN, VM_EVEN, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID),K_EVEN);
                                else
                                    SIGMA_EVEN = solve_transport_SIGMA_SC_3D(mat, EM_EVEN, VM_EVEN, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID),K_EVEN);
                                    SIGMA_EVEN.j = SIGMA_EVEN.j';
                                end

                                mat.n_of_modes=ceil(length(vecp1));
                                SIGMA_ODD = solve_transport_SIGMA_FC_3D(mat, EM_ODD, VM_ODD, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID),K_ODD);                                
                              
                                mat.n_of_modes=ceil(length(vecp1)+length(vecp2));
       

                                % mat.W_EVEN = SIGMA_EVEN.W;
                                mat.W_ODD  = SIGMA_ODD.W;
                                mat.toggleINIT = 1;
                                    
                                if mat.plot
                                    plotW(SIGMA_ODD.W, 13,2)
                                    % plotW(SIGMA_EVEN.W,24,2)

                                    plotW(SIGMA_ODD.RHO, 513,2)
                                    % plotW(SIGMA_EVEN.RHO,524,2)
                                end


                                n = SIGMA_ODD.n + SIGMA_EVEN.n;
                                
                                j = SIGMA_ODD.j + SIGMA_EVEN.j;
                            
                            elseif mat.doubleNm == 0
                                SIGMA = solve_transport_SIGMA_FC_3D(mat, EM, VM, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID),K);
                                n = SIGMA.n;
                                j = SIGMA.j;
                            end

                        end

                case 'DG'
                    tic
                    %% DG Method for solving in transport direction
                    DG = solve_transport_DG_3D(mat, EM, VM, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID));
                    toc
                    n = DG.n;
                    j = DG.j;
                    figure(999), plot(DG.chi, j);
                    results.Esub = EM;
                    results.Vsub = VM;

                    %% Map numeric discretization onto spatial discretization for Poisson
                    
                    n_sp = zeros(mat.Nx, mat.Ny, mat.Nz);
                    j_sp = zeros(mat.Nx);                   %yet not necessary

                    MapM = [0.5,0,0;0,1,0;0,0,0.5];

                    Mapping = zeros(mat.Nx,mat.Np);
                    for IN = 1 : N_chi
                        Mapping((IN-1)*(N_K_chi-1)+1:IN*(N_K_chi-1)+1 , (IN-1)*N_K_chi+1:IN*N_K_chi) = MapM;
                    end

                    Mapping(1,1)           = 1;
                    Mapping(mat.Nx,mat.Np) = 1;


                    for IY = 1 : mat.Ny
                        for IZ = 1 : mat.Nz
                            n_sp(:,IY,IZ) = Mapping*squeeze(n(:,IY,IZ));
                        end
                    end
                    results.chi = DG.chi;
                    results.xi  = DG.xi;
                    results.rho = DG.RHO;
                    results.p   = DG.p;
                    if (mat.dg.params.st_init == true)
                        results.rhs = DG.rhs;
                        results.rho = DG.rho;
                        results.A   = DG.A;
                    end
                    %%
                   



                    % jm = SIGMA.jm;
                % case 'OMICRON'
                % 
                %     OMICRON = solve_transport_OMICRON_2(mat, EM, VM, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID));
                % 
                %     n = OMICRON.nT;
                %     j = OMICRON.j;
                % 
                % case 'SIGMA_SCATTER'
                %     if performance_switch == 0
                %         [SIGMA, C_wo_drift] = solve_transport_S_SIGMA(mat, EM, VM, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID));
                %         performance_switch=1;
                %         results.C_wo_drift = C_wo_drift;
                %     else
                %         SIGMA = solve_transport_S_SIGMA_accel(mat, EM, VM, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID), C_wo_drift);
                %     end
                %     n = SIGMA.n;
                %     j = SIGMA.j;
                % 
                %     results.nm = SIGMA.nm;
                %     results.sigT = SIGMA.sigT;
                otherwise
                    
                    fprintf('no valid solver has been chosen\n consider: \n NEGF \t = Nonequilibriums Green Function Approach\n WIGNER\t = Wigner Transport Equation solver\n SIGMA \t = Liouville von Neumann Equation\n DG \t = Discontinous Galerkin Approach\n');
                    error('ending programm');
            end
            
            %disp(j(end-3));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % POISSON UPDATE

            if (strcmpi(solver,'DG')==true)
                [Vbi, iter_err] = solve_Poisson_3D(M, node_type, node_num, n_sp, mat, Vbi, Vg_bias);
            else
                [Vbi, iter_err] = solve_Poisson_3D(M, node_type, node_num, n, mat, Vbi, Vg_bias);
            end


            % if mat.plot
            %     figure(1),hold on, imagesc(squeeze(n(:,ceil(mat.Ny/2),:)));
            %     figure(2),hold on, imagesc(squeeze(Vbi(:,ceil(mat.Ny/2),:)));
            %     drawnow;
            % end


            iter_cnt = iter_cnt+1;
%             fprintf(fileID, '%.2f V at Drain running at %d iteration with accuracy %.2e\n', mat.Vd(ID), iter_cnt-1, iter_err);
            fprintf('%.2f V at Drain running at %d iteration with accuracy %.2e\n', mat.Vd(ID), iter_cnt-1, iter_err);
            fprintf('J end: %.6f J mean: %.6f with stdev: %.6f \n', j(end-3), mean(j(3:end-2)), std(j(3:end-2)));
        end       
        % end of self-consistent loop
        


        n_total(:,:,:, ID, IG) = n;
        j_total(:, ID, IG) = j;
        
        V_total(:,:,:, ID, IG) = Vbi;
        

    end
end



results.n_total = n_total;
results.j_total = j_total;
results.V_total = V_total;
results.Ef0=Ef0_L;

    

return


