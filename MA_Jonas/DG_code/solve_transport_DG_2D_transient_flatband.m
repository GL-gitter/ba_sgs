function DG = solve_transport_DG_2D_transient_flatband(mat, INIT)
    addpath DG_code/core/
    % clc
    % close all
    % hold on
    h = 6.626E-34;
    m0 = 9.109E-31;
    hb = h/2/pi;
    
    n_of_modes = mat.n_of_modes;
    n_of_valleys = mat.n_of_valleys;
    
    Esub_INIT = INIT.Esub;
    Vsub_INIT = INIT.Vsub;
    A_INIT    = INIT.A;   
    rhs_INIT  = INIT.rhs;
    rho_INIT  = INIT.rho;

    

    
    %% Get subband energies, system matrix and boundary conditions for gate voltage after switching 
    
    Vg = mat.Vg(2);
    
    DG = get_init_DG_2D_trans_flatband(mat, Vg);
    
    p    = DG.p;
    Esub = DG.Esub;
    Vsub = DG.Vsub;
    A    = DG.A;   
    rhs  = DG.rhs;
    rho  = DG.rho;
    
    %[A,rhs] = get_SysM(p);          % Get Stationary data
    %     rho = A\rhs;                    % Get Stationary solution
    % 
    %     p.U = 0.2;
    %     p.U1 = zeros(size(p.xg));
    %     p.U1(p.xg < p.L_chi/2+p.L_d/2 & p.xg > p.L_chi/2-p.L_d/2) = linspace(0, p.U, sum(p.xg < p.L_chi/2+p.L_d/2 & p.xg > p.L_chi/2-p.L_d/2));
    %     p.U1(p.xg >= p.L_chi/2+p.L_d/2) = p.U;
    % 
    %     [AU,rhsU] = get_SysM(p);        % Get Transient data
    % 
    %     if (p.doDG == true || p.doFV == true)
    %         rho_U = (reshape(rho, p.N_chi*p.N_K_chi, p.N_xi*p.N_K_xi)*p.phi.')';
    %         [n_over_t, j_over_t, rho_T] = transient(p, rhs, AU, rho, rho_U);
    %     elseif (p.doFEM == true)
    %         rho_U = (reshape(rho, p.N_chi*p.N_K_chi, p.N_xi*(p.N_K_xi-1)+1)*p.phi.')';
    %         [n_over_t, j_over_t, rho_T] = transient_FEM(p, rhs, AU, rho, rho_U);
    %     end
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
    
            %p = initParams(mat, EfL, EfR, EM, m, deg_factor(IV));
            %p.deg_factor_temp = deg_factor(IV);
    
            %% Mapping of Vsub, Esub onto DG_chi coordinates
            %Vsub_mod = zeros(p.N_chi*p.N_K_chi, mat.Ny);
            %m_mod    = zeros(p.N_chi*p.N_K_chi,1);
            %Esub_mod = zeros(p.N_chi*p.N_K_chi, mat.Ny, mat.Nz);
            
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
    
    %% stepping through time
    n_over_t = zeros(p.N_chi*p.N_K_chi, mat.Ny, mat.Nt);
    j_over_t = zeros(p.N_chi*p.N_K_chi, mat.Nt);
    
    n_over_t_mode = zeros(n_of_modes, p.N_chi*p.N_K_chi, mat.Ny, mat.Nt);
    j_over_t_mode = zeros(n_of_modes, p.N_chi*p.N_K_chi, mat.Nt);
    fid = figure('name', "Dichte");
    
    Nt = mat.Nt;
    dt = mat.dt;
    %% Try to depack A,rho and rhs for faster depacking
    % sadly no faster option here

    %% Transiente Berechnung mit Runge-Kutta 4. Ordnung
    for tstep = 1 : Nt
        if mod((tstep-1),5) == 0
            display(['Zeitschritt ', num2str(tstep), '/', num2str(Nt)])
        end
        parfor IM = 1:n_of_modes
            % Runge-Kutta 4. Ordnung
            tic
            k1 = rhs{IM}-A{IM}*rho_INIT{IM};
            k2 = rhs{IM}-A{IM}*(rho_INIT{IM}+dt/2*k1);
            k3 = rhs{IM}-A{IM}*(rho_INIT{IM}+dt/2*k2);
            k4 = rhs{IM}-A{IM}*(rho_INIT{IM}+dt*k3);

            rho_INIT{IM} = rho_INIT{IM}+(dt/6)*(k1+2*k2+2*k3+k4);

            result_matrix = reshape(squeeze(rho_INIT{IM}), p.N_chi*p.N_K_chi, p.Nk)*p.phi.';
            toc
           
            % Observablen speichern
            %kron(ones(1,p.N_xi), real((result_matrix(:, p.N_xi*p.N_K_xi/2+1)+ result_matrix(:, p.N_xi*p.N_K_xi/2))/2)).*reshape(abs(Vsub_mod(IM, :, :)).^2, p.N_chi*p.N_K_chi, mat.Ny)/mat.dy/1E-9;
            if mod(p.N_xi,2)==0
                n_over_t_mode(IM, :, :, tstep) = kron(ones(1,mat.Ny), real((result_matrix(:, (p.N_xi)/2+1)+ result_matrix(:, (p.N_xi)/2))/2)).*reshape(abs(Vsub_mod(IM, :, :)).^2, p.N_chi*p.N_K_chi, mat.Ny)/mat.dy/1E-9;   
                j_over_t_mode(IM, :, tstep) = p.q*hb./m_mod.*imag(result_matrix(:, (p.N_xi)/2+1)-result_matrix(:, (p.N_xi)/2))/p.delta_xi; 
            elseif mod(p.N_xi,2)==1
                n_over_t_mode(IM, :, :, tstep) = kron(ones(1,mat.Ny), real(result_matrix(:, (p.N_xi*p.N_K_xi-1)/2+1))).*reshape(abs(Vsub_mod(IM, :, :)).^2, p.N_chi*p.N_K_chi, mat.Ny)/mat.dy/1E-9;   
                j_over_t_mode(IM, :, tstep) = p.q*hb./m_mod.*imag(result_matrix(:, (p.N_xi-1)/2+2)-result_matrix(:, (p.N_xi-1)/2))/p.delta_xi; 
            end                           
        end

        n_over_t(:,:,tstep) = sum(squeeze(n_over_t_mode(:,:,:,tstep)), 1);
        %n_over_t(:,:,tstep) = squeeze(n_over_t_mode(1, :, :, tstep))+squeeze(n_over_t_mode(2, :, :, tstep))+squeeze(n_over_t_mode(3, :, :, tstep));
        j_over_t(:,tstep) = sum(squeeze(j_over_t_mode(:,:,tstep)), 1);

        if mod((tstep-1),5) == 0
            figure(fid)
            title(['t=', num2str(mat.t(tstep)), 'fs']);
            plot(p.chi_Koord,real(squeeze(n_t_pre(:,ceil(mat.Ny/2)))), 'r');
            hold on
            plot(p.chi_Koord,real(squeeze(n_over_t(:,ceil(mat.Ny/2),tstep))), 'b');
            %plot(p.chi_Koord,real(result_matrix(:, p.N_xi*p.N_K_xi/2)), 'r');
            hold off
        end
    end
    
    %%

     DG.n_over_t = n_over_t(:,:,1:10:end);
     DG.j_over_t = j_over_t(:,1:10:end);
     DG.rho_INIT = rho_INIT;
     DG.chi = p.chi_Koord;
     DG.xi  = p.xi_Rechengebiet_Elem;
end






  

