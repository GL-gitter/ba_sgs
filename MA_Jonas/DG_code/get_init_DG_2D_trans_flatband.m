function INIT = get_init_DG_2D_trans_flatband(mat, Vg)
    %addpath DG_code/core/
    %% Adaption for DG algorithm
    N_K_chi = mat.dg.params.N_K_chi;
    N_chi   = ((mat.Nx-1)/(N_K_chi-1));
    mat.Np  = N_chi*N_K_chi;
    %%
    
    %solve for the subbands
    [Esub, ~] = solve_subbands(mat);
    
    %solve for the fermi level
    [Ef0_L, Ef0_R] = solve_fermi(Esub, mat);

    %eval poisson
    %[M, node_num, node_type] = eval_Poisson(mat);   %not needed for flatband transient
    
    Vg_bias = Ef0_L + mat.phi_m_g-mat.Xi-Vg;
    
    %set Vbi
    Vbi = Vg_bias*repmat(mat.boundary(:,1), 1, mat.Ny);
    Vbi(mat.x>mat.L_s+mat.L_c,:) = -mat.Vd(1);
    Vbi=sgolayfilt(Vbi,0,17);
    Vbi(((mat.x>mat.L_s) .* (mat.x<mat.L_s+mat.L_c))==1,1)   = Vg_bias;
    Vbi(((mat.x>mat.L_s) .* (mat.x<mat.L_s+mat.L_c))==1,end) = Vg_bias;

    Vref = mat.V; 
    mat.V = Vref+Vbi;

    [EM, VM] = solve_subbands(mat);
    
    for IM = 1:mat.n_of_modes
        EM(1, :, IM)=squeeze(sgolayfilt(EM(1,:,IM),0,mat.smoothEM*2+1));
        VM(1,:,IM,:)=sgolayfilt(squeeze(VM(1,:,IM,:)),0,mat.smoothVM*2+1);

        for IX=1:mat.Nx
            VM(1,IX,IM,:) = real(VM(1,IX,IM,:))/sqrt(sum(real(VM(1,IX,IM,:)).^2,'all')); %need to renormalize WFs
        end
    end

    mat.V = Vref;   %not quiet sure if needed right here
   
    %% INTERFACE
    %solve_transport_DG_2D(mat, EM, VM, Ef0_L-mat.Vs, Ef0_R-mat.Vd(ID));
    %solve_transport_DG_2D_transient(mat, Esub, Vsub, EfL, EfR)  function head
    Esub = EM;
    Vsub = VM;
    EfL  = Ef0_L-mat.Vs;
    EfR  = Ef0_R-mat.Vd(1);
    %%
    
    h = 6.626E-34;
    m0 = 9.109E-31;
    hb = h/2/pi;
    
    n_of_valleys = mat.n_of_valleys;
    n_of_modes = mat.n_of_modes;
    deg_factor = mat.deg_factor;
    
    
    N_K_chi = mat.dg.params.N_K_chi;
    
    
    %% initParams structure needed for parfor loop (make more efficient in future)
    if (n_of_modes==1)
        p     = initParams(mat, EfL, EfR, squeeze(Esub(1, :, 1)), squeeze(mat.me_x(1, :, ceil(mat.Ny/2)))*m0, deg_factor(1));
    else
        p     = initParams(mat, EfL, EfR, squeeze(Esub(1, :, 1)), squeeze(mat.me_x(1, :, ceil(mat.Ny/2)))*m0, deg_factor(1));
        p_top = repmat(p,n_of_valleys,n_of_modes); 
    
        for IV=1:n_of_valleys
            for IM=1:n_of_modes                      %% hast to be reworked for moren than on valley (n_of_valley > 1)!
                p = initParams(mat, EfL, EfR, squeeze(Esub(IV, :, IM)), squeeze(mat.me_x(IV, :, ceil(mat.Ny/2)))*m0, deg_factor(IV));  
                p_top(IV,IM) = p;
            end
        end
    end
    % A_col   = struct;
    % rhs_col = struct;
    % rho_col = struct;
    Atemp = cell(1,n_of_modes);
    rhstemp = cell(1,n_of_modes);
    rhotemp = cell(1,n_of_modes);
    
    for IV=1:n_of_valleys
        for IM=1:n_of_modes
            
            %% Valley and Mode specific values
            %EM = squeeze(Esub(IV, :, IM));                           % Valley and Mode specific Subband energy
            m = squeeze(mat.me_x(IV, :, ceil(mat.Ny/2)))*m0;
            %%
            if (n_of_modes>1)
              p = p_top(IV,IM);
            end
            %p = initParams(mat, EfL, EfR, EM, m, deg_factor(IV));
            %p.deg_factor_temp = deg_factor(IV);

            %% Mapping of Vsub, Esub onto DG_chi coordinates
            Vsub_mod = zeros(p.N_chi*p.N_K_chi, mat.Ny);
            m_mod    = zeros(p.N_chi*p.N_K_chi,1);
            %Esub_mod = zeros(p.N_chi*p.N_K_chi, mat.Ny, mat.Nz);
            for Nchi=1:p.N_chi
                V_temp = reshape(Vsub(IV,(Nchi-1)*(p.N_K_chi-1)+1:(Nchi)*(p.N_K_chi-1)+1, IM,:),p.N_K_chi,mat.Ny);     
                Vsub_mod((Nchi-1)*p.N_K_chi+1:Nchi*p.N_K_chi,:) = V_temp;
    
                m_temp = reshape(m((Nchi-1)*(p.N_K_chi-1)+1:(Nchi)*(p.N_K_chi-1)+1),p.N_K_chi,1);
                m_mod((Nchi-1)*p.N_K_chi+1:Nchi*p.N_K_chi) = m_temp;
            end
            
            %m_mod = m_mod.';
            %%
            [A,rhs] = get_SysM(mat, p, squeeze(Esub(IV, :, IM)));
            rhs = sparse(rhs);
            %% Try scaling factor for A and rhs!!! 
            rhs =(1e-11)*rhs;      
            A   =(1e-11)*A;
            rho = A\rhs;

            Atemp{IM} = A;
            rhstemp{IM} = rhs;
            rhotemp{IM} = rho;
        end
    end


INIT.A   = Atemp;
INIT.rhs = rhstemp;
INIT.rho = rhotemp;

INIT.Esub = Esub;
INIT.Vsub = Vsub;
INIT.p    = p;

return