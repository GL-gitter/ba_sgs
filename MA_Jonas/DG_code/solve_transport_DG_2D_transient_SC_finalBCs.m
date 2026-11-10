function DG = solve_transport_DG_2D_transient_SC_finalBCs(mat, INIT)
    addpath DG_code/core/
    % clc
    % close all
    % hold on
    h = 6.626E-34;
    m0 = 9.109E-31;
    hb = h/2/pi;

    N_K_chi = mat.dg.params.N_K_chi;
    N_chi   = ((mat.Nx-1)/(N_K_chi-1));
    mat.Np      = N_chi*N_K_chi;

    deg_factor = mat.deg_factor;
    n_of_modes = mat.n_of_modes;
    n_of_valleys = mat.n_of_valleys;
    
    Esub_INIT = INIT.Esub;
    Vsub_INIT = INIT.Vsub;
    A_INIT    = INIT.A;   
    rhs_INIT  = INIT.rhs;
    rho_INIT  = INIT.rho;
    p         = INIT.p;

    Vbi = INIT.V_total;
  
    %% Get subband energies, system matrix and boundary conditions for gate voltage after switching 
       
    Vg_temp = mat.Vg;
    mat.Vg  = Vg_temp(end);
    DG = solve_self_consistent(mat, 'DG');

    mat.Vg = Vg_temp;
    % 
    % p    = DG.p;
    % Esub = DG.Esub;
    % Vsub = DG.Vsub;
    % A    = DG.A;
    rhs = cell(1,n_of_modes);
    for IM = 1:n_of_modes
     rhs{IM}  = mat.dg.params.MPC*DG.rhs{IM};
    end
    % rho  = DG.rho;
    
  
    rho_pre = zeros(n_of_modes, p.N_chi*p.N_K_chi, p.N_xi);

    for IM = 1:n_of_modes
        rho_pre(IM,:,:) = (reshape(squeeze(rho_INIT{IM}), p.N_chi*p.N_K_chi, p.Nk)*p.phi.');
    end
    
    %% preparing eigenvectors of modes for calculation of carrier density
    Vsub_mod     = zeros(n_of_modes, p.N_chi*p.N_K_chi, mat.Ny);
    Vsub_mod_pre = zeros(n_of_modes, p.N_chi*p.N_K_chi, mat.Ny);
    m_mod        = zeros(n_of_modes, p.N_chi*p.N_K_chi,1);
    n_pre        = zeros(n_of_modes, p.N_chi*p.N_K_chi, mat.Ny);
    
    for IM=1:n_of_modes            
            %% Valley and Mode specific values
            %EM = squeeze(Esub(IV, :, IM));                           % Valley and Mode specific Subband energy
            m = squeeze(mat.me_x(1, :, ceil(mat.Ny/2)))*m0;
   
            %% Mapping of Vsub, Esub onto DG_chi coordinates            
            for Nchi=1:p.N_chi
                V_temp = reshape(Vsub_INIT(1,(Nchi-1)*(p.N_K_chi-1)+1:(Nchi)*(p.N_K_chi-1)+1, IM,:),p.N_K_chi,mat.Ny);     
                Vsub_mod(IM, (Nchi-1)*p.N_K_chi+1:Nchi*p.N_K_chi,:) = V_temp;

                V_temp_pre = reshape(Vsub_INIT(1,(Nchi-1)*(p.N_K_chi-1)+1:(Nchi)*(p.N_K_chi-1)+1, IM,:),p.N_K_chi,mat.Ny);     
                Vsub_mod_pre(IM, (Nchi-1)*p.N_K_chi+1:Nchi*p.N_K_chi,:) = V_temp_pre;
    
                m_temp = reshape(m((Nchi-1)*(p.N_K_chi-1)+1:(Nchi)*(p.N_K_chi-1)+1),p.N_K_chi,1);
                m_mod(IM, (Nchi-1)*p.N_K_chi+1:Nchi*p.N_K_chi) = m_temp;
            end
            
            if mod(p.N_xi,2)==0
                n_pre(IM, :, :) = kron(ones(1,mat.Ny), (squeeze(real(squeeze(rho_pre(IM, :, p.N_xi*p.N_K_xi/2+1))+ squeeze(rho_pre(IM, :, p.N_xi*p.N_K_xi/2))))/2).').*reshape(abs(Vsub_mod_pre(IM, :, :)).^2, p.N_chi*p.N_K_chi, mat.Ny)/mat.dy/1E-9;   
            elseif mod(p.N_xi,2)==1
                n_pre(IM, :, :) = kron(ones(1,mat.Ny), squeeze(real(rho_pre(IM,:, (p.N_xi*p.N_K_xi-1)/2+1))).').*reshape(abs(Vsub_mod_pre(IM, :, :)).^2, p.N_chi*p.N_K_chi, mat.Ny)/mat.dy/1E-9;   
            end     
    end
    
    n_t_pre = squeeze(sum(n_pre,1));

    m_mod = squeeze(m_mod(1,:,:));
    m_mod = m_mod.';
    


    %% Stepping through time
    n_over_t = zeros(p.N_chi*p.N_K_chi, mat.Ny, mat.Nt);
    j_over_t = zeros(p.N_chi*p.N_K_chi, mat.Nt);
    
    n_over_t_mode = zeros(n_of_modes, p.N_chi*p.N_K_chi, mat.Ny, mat.Nt);
    j_over_t_mode = zeros(n_of_modes, p.N_chi*p.N_K_chi, mat.Nt);
    fid = figure('name', "Dichte");
    
    Nt = mat.Nt;
    dt = mat.dt;
    %% Try to depack A,rho and rhs for faster depacking
    % sadly no faster option here
    %% Prepare mapping matrix to map DG discetization of carrier density onto FV discretization for Poisson
    n_sp = zeros(mat.Nx, mat.Ny);
    %j_sp = zeros(mat.Nx);                   %yet not necessary

    MapM = [0.5,0,0;0,1,0;0,0,0.5];

    Mapping = zeros(mat.Nx,mat.Np);
    for IN = 1 : N_chi
        Mapping((IN-1)*(N_K_chi-1)+1:IN*(N_K_chi-1)+1 , (IN-1)*N_K_chi+1:IN*N_K_chi) = MapM;
    end

    Mapping(1,1)           = 1;
    Mapping(mat.Nx,mat.Np) = 1;

    IV = 1;

    A = cell(1,n_of_modes);
    rhs = cell(1,n_of_modes);
    %% Transiente Berechnung mit Runge-Kutta 4. Ordnung
    for tstep = 1 : Nt
        tic
        if mod((tstep-1),5) == 0
            display(['Zeitschritt ', num2str(tstep), '/', num2str(Nt)])
        end
        
        %% Refresh potential for calculation of new subband energy
        Vref = mat.V; mat.V = Vref+Vbi;
        [EM, VM] = solve_subbands(mat);
        mat.V = Vref;

        [Ef0_L, Ef0_R] = solve_fermi(EM, mat);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % Unsure about boundary conditions, recalculate in every timestep
        % or approximation of time independent BCs sufficient 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Runge-Kutta 4. Ordnung
        for IM = 1:n_of_modes
            p = initParams(mat, Ef0_L-mat.Vs, Ef0_R-mat.Vd, squeeze(EM(IV, :, IM)), m, deg_factor(IV));   %% Mods needed when simulating time dependent Vd; Mods needed if n_of_valleys > 1
            [A{IM},~] = get_SysM(mat, p, squeeze(EM(IV, :, IM)));
            
            A{IM}   = (2.7e-12)*A{IM};
            %rhs{IM} = (3e-12)*rhs{IM};
            
            k1 = rhs{IM}-A{IM}*rho_INIT{IM};
            k2 = rhs{IM}-A{IM}*(rho_INIT{IM}+dt/2*k1);
            k3 = rhs{IM}-A{IM}*(rho_INIT{IM}+dt/2*k2);
            k4 = rhs{IM}-A{IM}*(rho_INIT{IM}+dt*k3);

            rho_INIT{IM} = rho_INIT{IM}+(dt/6)*(k1+2*k2+2*k3+k4);

            result_matrix = reshape(squeeze(rho_INIT{IM}), p.N_chi*p.N_K_chi, p.Nk)*p.phi.';
        
            
            % Store charge carrier density and current density
            %kron(ones(1,p.N_xi), real((result_matrix(:, p.N_xi*p.N_K_xi/2+1)+ result_matrix(:, p.N_xi*p.N_K_xi/2))/2)).*reshape(abs(Vsub_mod(IM, :, :)).^2, p.N_chi*p.N_K_chi, mat.Ny)/mat.dy/1E-9;
            if mod(p.N_xi,2)==0
                n_over_t_mode(IM, :, :, tstep) = kron(ones(1,mat.Ny), real((result_matrix(:, (p.N_xi)/2+1)+ result_matrix(:, (p.N_xi)/2))/2)).*reshape(abs(Vsub_mod(IM, :, :)).^2, p.N_chi*p.N_K_chi, mat.Ny)/mat.dy/1E-9;   
                j_over_t_mode(IM, :, tstep) = p.q*hb./m_mod.*imag(result_matrix(:, (p.N_xi)/2+1)-result_matrix(:, (p.N_xi)/2))/p.delta_xi; 
            elseif mod(p.N_xi,2)==1
                n_over_t_mode(IM, :, :, tstep) = kron(ones(1,mat.Ny), real(result_matrix(:, (p.N_xi*p.N_K_xi-1)/2+1))).*reshape(abs(Vsub_mod(IM, :, :)).^2, p.N_chi*p.N_K_chi, mat.Ny)/mat.dy/1E-9;   
                j_over_t_mode(IM, :, tstep) = p.q*hb./m_mod.*imag(result_matrix(:, (p.N_xi-1)/2+2)-result_matrix(:, (p.N_xi-1)/2))/p.delta_xi; 
            end
                
        end
        if (n_of_modes > 1)
            n_over_t(:,:,tstep) = sum(squeeze(n_over_t_mode(:,:,:,tstep)), 1);
            %n_over_t(:,:,tstep) = squeeze(n_over_t_mode(1, :, :, tstep))+squeeze(n_over_t_mode(2, :, :, tstep))+squeeze(n_over_t_mode(3, :, :, tstep));
            j_over_t(:,tstep) = sum(squeeze(j_over_t_mode(:,:,tstep)), 1);
        else
            n_over_t(:,:,tstep) = squeeze(n_over_t_mode(:,:,:,tstep));
            %n_over_t(:,:,tstep) = squeeze(n_over_t_mode(1, :, :, tstep))+squeeze(n_over_t_mode(2, :, :, tstep))+squeeze(n_over_t_mode(3, :, :, tstep));
            j_over_t(:,tstep) = squeeze(j_over_t_mode(:,:,tstep));
        end

        if mod((tstep-1),5) == 0
            figure(fid)
            title(['t=', num2str(mat.t(tstep)), 'fs']);
            plot(p.chi_Koord,real(squeeze(n_t_pre(:,ceil(mat.Ny/2)))), 'r');
            hold on
            plot(p.chi_Koord,real(squeeze(n_over_t(:,ceil(mat.Ny/2),tstep))), 'b');
            %plot(p.chi_Koord,real(result_matrix(:, p.N_xi*p.N_K_xi/2)), 'r');
            hold off
        end
        % solve poisson
        % refresh Subbands
        % build new Systemmatrix
        % calculate new time step
        for IY = 1 : mat.Ny
             n_sp(:,IY) = Mapping*squeeze(n_over_t(:,IY,tstep));
        end

        Vbi_IN = INIT.V_total;
        Vbi = solve_Poisson2D_transient(mat, Vbi_IN, INIT.Ef0, tstep, n_sp);
        toc
    end
    
    %%

     DG.n_over_t = n_over_t(:,:,1:10:end);
     DG.j_over_t = j_over_t(:,1:10:end);
     DG.rho_INIT = rho_INIT;
     DG.chi = p.chi_Koord;
     DG.xi  = p.xi_Rechengebiet_Elem;
end






  

