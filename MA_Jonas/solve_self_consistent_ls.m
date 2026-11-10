function results = solve_self_consistent_ls(mat, solver)
% For scattering please use non load/save version!

if mat.firstrun
    mat.time_start = datetime;
    % solve for the subbands and fermi
    [Esub, ~] = solve_subbands(mat);
    [Ef0_L, Ef0_R] = solve_fermi(Esub, mat);
    % init matrices
    [M, node_num, node_type] = eval_Poisson(mat);
    n_total = zeros(mat.Nx, mat.Ny, length(mat.Vd));
    j_total = zeros(mat.Nx, length(mat.Vd));
    V_total = zeros(mat.Nx, mat.Ny, length(mat.Vd));
    performance_switch = 0;
    
    IG0=1;
    ID0=1;
            
    iter_cnt = 0;
    iter_err = 1;
    
else
    load(strcat('/work/smmapech/msa_matlab/inits/', string(mat.job_id),'_mat.mat') );
    mat.time_start = datetime;
    %mat.time_max   = 2; % in minutes
    IG0=IG;
    ID0=ID;
    aaaa=1;
    
end
        

for IG = IG0 : length(mat.Vg)
    
    
    Vg_bias = Ef0_L + mat.phi_m_g-mat.Xi-mat.Vg(IG); % -mat.Vg(IG); % Matched work functions  v + mat.phi_m_g-mat.Xi-mat.Vg(IG);
    
    CPL=0;
    for ID = ID0 : length(mat.Vd)
        if iter_cnt==0
            if ID > 2
                
                Vbi = 2*V_total(:,:,ID-1)-V_total(:,:,ID-2);
                
            else
                
                Vbi = Vg_bias*repmat(mat.boundary(:,1), 1, mat.Ny);
                
                Vbi(mat.x>mat.L_s+mat.L_c,:) = -mat.Vd(ID);
                Vbi=sgolayfilt(Vbi,0,17);
                Vbi(((mat.x>mat.L_s) .* (mat.x<mat.L_s+mat.L_c))==1,1)   =Vg_bias;
                Vbi(((mat.x>mat.L_s) .* (mat.x<mat.L_s+mat.L_c))==1,end) =Vg_bias;
            end
        end
        
        
        while iter_err > mat.poisson.opt.iter_err && iter_cnt < mat.poisson.opt.iter_max
            

            abs(minutes(mat.time_start-datetime))
            if abs(minutes(mat.time_start-datetime)) > mat.time_max
                    fprintf('Saving temporary results')
                    save(strcat('/work/smmapech/msa_matlab/inits/', string(mat.job_id),'_mat.mat') );
                    fprintf('Done, return to main');
                    results.finished = false;
                    return
            end


            % if length(mat.Vd)==1
            % 
            %     if iter_cnt>0 & mat.CPL>0 & ~strcmp(solver, 'NEGF') & ~strcmp(solver, 'NEGF_RSMV')
            %         CPL=1;
            %         disp('Self coupling turned on');
            %     end
            % 
            %     if iter_cnt>10 & mat.CPL==2 & ~strcmp(solver, 'NEGF') & ~strcmp(solver, 'NEGF_RSMV')
            %         CPL=2;
            %         disp('Full coupling turned on');
            %     end
            % 
            % else
            %     if ID < 2 & mat.CPL>0 & ~strcmp(solver, 'NEGF') & ~strcmp(solver, 'NEGF_RSMV')
            %         CPL=1;
            %         disp('Self coupling turned on');
            %     elseif ID>2 & mat.CPL>1 & ~strcmp(solver, 'NEGF') & ~strcmp(solver, 'NEGF_RSMV')  & iter_cnt>20
            %         CPL=2;
            %         disp('Full coupling turned on');
            %     end
            % end

            Vref = mat.V; mat.V = Vref+Vbi;
            [EM, VM] = solve_subbands(mat);
            
            for IM = 1:mat.n_of_modes
                    EM(1, :, IM)=squeeze(sgolayfilt(EM(1,:,IM),0,mat.smoothEM*2+1));
                    VM(1,:,IM,:)=sgolayfilt(squeeze(VM(1,:,IM,:)),0,mat.smoothVM*2+1);

                    for IX=1:mat.Nx
                        VM(1,IX,IM,:) = real(VM(1,IX,IM,:))/sqrt(sum(real(VM(1,IX,IM,:)).^2,'all')); %need to renormalize WFs
                    end
                end
            
            
            
            
            K = getK_alt(mat.Nx,mat.dx*1e-9,mat.n_of_modes, 1,mat.me_x_ch, mat.me_x, VM, mat.K_order);
            
            if mat.plot; plotK(K,44,mat.n_of_modes); end
            
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
                
                case 'NEGF'
                    NEGF = solve_transport_NEGF(mat, EM, VM, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID));
                    
                    n = NEGF.n;
                    j = NEGF.j;

                % case 'NEGF_RS' % discontinued
                %         NEGF = solve_transport_NEGF_real_space(mat, Vref+Vbi, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID));
                %     end
    
                case 'NEGF_RSMV'
                    tic
                    NEGF_RSMV = solve_transport_NEGF_real_space_MV(mat, Vref+Vbi, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID));
                    toc
                    n = NEGF_RSMV.n;
                    j = NEGF_RSMV.j;
                    results.n2=NEGF_RSMV.n2;
                    results.j2=NEGF_RSMV.j2;
                
                case 'NEGF_CMS'
                    tic
                    NEGF_CMS = solve_transport_NEGF_CMS_MV(mat, Vref+Vbi, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID),VM);
                    toc
                        
                    n = NEGF_CMS.n;
                    j = NEGF_CMS.j;
                    results.n2=NEGF_CMS.n2;
                    results.j2=NEGF_CMS.j2;

                case 'NEGF_CMS_ALT'
                    
                    tic
                    NEGF_CMS_ALT = solve_transport_NEGF_CMS_alt(mat, Vref+Vbi, EM, VM, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID),K);
                    toc

                    n = NEGF_CMS_ALT.n;
                    j = NEGF_CMS_ALT.j;
                    results.n2=NEGF_CMS_ALT.n2;
                    results.j2=NEGF_CMS_ALT.j2;

                case 'WIGNER'
                    
                    if performance_switch == 0
                        
     
                        WIGNER = solve_transport_WIGNERconst_eff_optim(mat, EM, VM, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID));
                        results.C_wo_drift = WIGNER.C_wo_drift;
                        performance_switch = 1;


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
                            %SIGMA = solve_transport_SIGMA(mat, EM, VM, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID));
                            SIGMA = solve_transport_SIGMA_SC(mat, EM, VM, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID),K);
                            n = SIGMA.n;
                            j = SIGMA.j;
                            results.wf = SIGMA.wf;
                        
                        elseif CPL == 1
                            SIGMA = solve_transport_SIGMA_SC(mat, EM, VM, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID),K);
                            n = SIGMA.n;
                            j = SIGMA.j;
                            results.wf = SIGMA.wf;

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
        
                                VM_ODD  = VM(1,:,vecp1,:);
                                VM_EVEN = VM(1,:,vecp2,:);

                                
                                K_ODD.K0 = K.K0(vecp1,vecp1,:);
                                K_ODD.K1 = K.K1(vecp1,vecp1,:);
                                K_ODD.K2 = K.K2(vecp1,vecp1 ,:);

                                K_EVEN.K0 = K.K0(vecp2,vecp2,:);
                                K_EVEN.K1 = K.K1(vecp2,vecp2,:);
                                K_EVEN.K2 = K.K2(vecp2,vecp2,:);

                                mat.n_of_modes=ceil(length(vecp2));

                                if mat.n_of_modes>1
                                    SIGMA_EVEN = solve_transport_SIGMA_FC(mat, EM_EVEN, VM_EVEN, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID),K_EVEN);
                                else
                                    SIGMA_EVEN = solve_transport_SIGMA_SC(mat, EM_EVEN, VM_EVEN, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID),K_EVEN);
                                    SIGMA_EVEN.j = SIGMA_EVEN.j';
                                end

                                mat.n_of_modes=ceil(length(vecp1));
                                SIGMA_ODD = solve_transport_SIGMA_FC(mat, EM_ODD, VM_ODD, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID),K_ODD);                                
                              
                                mat.n_of_modes=ceil(length(vecp1)+length(vecp2));
       

                                % mat.W_EVEN = SIGMA_EVEN.W;
                                % mat.W_ODD  = SIGMA_ODD.W;
                                mat.toggleINIT = 1;
                                    
                                if mat.plot
                                    plotW(SIGMA_ODD.W, 13,2)
                                    % plotW(SIGMA_EVEN.W,24,2)

                                    plotW(SIGMA_ODD.RHO, 513,2)
                                    % plotW(SIGMA_EVEN.RHO,524,2)
                                end


                                n = SIGMA_ODD.n + SIGMA_EVEN.n;
                                
                                j = SIGMA_ODD.j + SIGMA_EVEN.j;

                                results.W_ODD=SIGMA_ODD.W;
                                results.W_EVEN=SIGMA_EVEN.W;

                            elseif mat.doubleNm == 0
                                SIGMA = solve_transport_SIGMA_FC(mat, EM, VM, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID),K);
                                n = SIGMA.n;
                                j = SIGMA.j;
                                results.W=SIGMA.W;
                            end

                        end
                    
                    
                    
                    % jm = SIGMA.jm;
                case 'OMICRON'
                    
                    OMICRON = solve_transport_OMICRON_2(mat, EM, VM, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID));
                    
                    n = OMICRON.nT;
                    j = OMICRON.j;

                case 'SIGMA_SCATTER'
                    if performance_switch == 0
                        [SIGMA, C_wo_drift] = solve_transport_S_SIGMA(mat, EM, VM, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID));
                        performance_switch=1;
                        results.C_wo_drift = C_wo_drift;
                    else
                        SIGMA = solve_transport_S_SIGMA_accel(mat, EM, VM, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID), C_wo_drift);
                    end
                    n = SIGMA.n;
                    j = SIGMA.j;
                    
                    results.nm = SIGMA.nm;
                    results.sigT = SIGMA.sigT;
                otherwise
                    
                    fprintf('no valid solver has been chosen\n consider: \n NEGF \t = Nonequilibriums Green Function Approach\n WIGNER\t = Wigner Transport Equation solver\n SIGMA \t = Liouville von Neumann Equation\n');
                    error('ending programm');
            end
            %disp(j(end-3));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % POISSON UPDATE
            if mat.plot
                figure(2),hold on, plot(squeeze(n(:,ceil(mat.Ny/2))));
                
                figure(1),mesh(n)
                drawnow;
            end
            [Vbi, iter_err] = solve_Poisson(M, node_type, node_num, n, mat, Vbi, Vg_bias);
            

            iter_cnt = iter_cnt+1;
%             fprintf(fileID, '%.2f V at Drain running at %d iteration with accuracy %.2e\n', mat.Vd(ID), iter_cnt-1, iter_err);
            fprintf('%.2f V at Drain running at %d iteration with accuracy %.2e\n', mat.Vd(ID), iter_cnt-1, iter_err);
            fprintf('J end: %.6f J mean: %.6f with stdev: %.6f \n', j(end-3), mean(j(3:end-2)), std(j(3:end-2)));
        end

        iter_cnt = 0;
        iter_err = 1;

        n_total(:,:, ID) = n;
        j_total(:, ID) = j;
        V_total(:,:, ID) = Vbi;        

    end
    
end

results.n_total = n_total;
results.j_total = j_total;
results.V_total = V_total;
results.finished = true;
results.Ef0=Ef0_L;
