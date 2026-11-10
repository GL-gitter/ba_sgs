function DG = solve_transport_DG_3D(mat, Esub, Vsub, EfL, EfR)

% clc
% close all
% hold on
h = 6.626E-34;
m0 = 9.109E-31;
hb = h/2/pi;

addpath DG_code/core/

DG = struct;

n_of_valleys = mat.n_of_valleys;
n_of_modes = mat.n_of_modes;
deg_factor = mat.deg_factor;


N_K_chi = mat.dg.params.N_K_chi;
nT      = zeros(((mat.Nx-1)/(N_K_chi-1))*N_K_chi, mat.Ny, mat.Nz);
jT      = zeros(((mat.Nx-1)/(N_K_chi-1))*N_K_chi,1);
RHO_DG  = zeros(n_of_valleys,n_of_modes,((mat.Nx-1)/(N_K_chi-1))*N_K_chi,mat.dg.params.N_xi);

n_mode  = zeros(n_of_valleys, n_of_modes, ((mat.Nx-1)/(N_K_chi-1))*N_K_chi, mat.Ny, mat.Nz);
j_mode  = zeros(n_of_valleys, n_of_modes, ((mat.Nx-1)/(N_K_chi-1))*N_K_chi,1);
IN_max  = size(n_mode,3);   %assignment required for parfor loop


%% initParams structure needed for parfor loop (make more efficient in future)
p     = initParams(mat, EfL, EfR, squeeze(Esub(1, :, 1)), squeeze(mat.me_x(1, :, ceil(mat.Ny/2),ceil(mat.Nz/2)))*m0, deg_factor(1));
p_top = repmat(p,n_of_valleys,n_of_modes); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mat.MP ==1
    m_temp = repmat(mat.mexv,1,mat.Nx);
else
    m_temp = repmat(squeeze(mat.me_x(1, :, ceil(mat.Ny/2),ceil(mat.Nz/2))),n_of_modes,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for IV=1:n_of_valleys
    for IM=1:n_of_modes                      %% hast to be reworked for moren than on valley (n_of_valley > 1)!
        p = initParams(mat, EfL, EfR, squeeze(Esub(IV, :, IM)), squeeze(m_temp(IM,:))*m0, deg_factor(IV));
        p_top(IV,IM) = p;
    end
end

if (mat.dg.params.st_init == true)
    Atemp = cell(1,n_of_modes);
    rhstemp = cell(1,n_of_modes);
    rhotemp = cell(1,n_of_modes); 

end

for IV=1:n_of_valleys
    for IM=1:n_of_modes
%%

        %% Valley and Mode specific values
        %EM = squeeze(Esub(IV, :, IM));    
        % Valley and Mode specific Subband energy
        if mat.MP ==1
            m_temp = repmat(mat.mexv,1,mat.Nx);
        else
            m_temp = repmat(squeeze(mat.me_x(1, :, ceil(mat.Ny/2),ceil(mat.Nz/2))),n_of_modes,1);
        end
        m = squeeze(m_temp(IM,:))*m0;
        
        %%
        p = p_top(IV,IM);
        %p = initParams(mat, EfL, EfR, EM, m, deg_factor(IV));
        %p.deg_factor_temp = deg_factor(IV);

        %% Mapping of Vsub, Esub onto DG_chi coordinates
        Vsub_mod = zeros(p.N_chi*p.N_K_chi, mat.Ny, mat.Nz);
        m_mod    = zeros(p.N_chi*p.N_K_chi,1);
        %Esub_mod = zeros(p.N_chi*p.N_K_chi, mat.Ny, mat.Nz);
        
        for Nchi=1:p.N_chi
            V_temp = reshape(Vsub(IV,(Nchi-1)*(p.N_K_chi-1)+1:(Nchi)*(p.N_K_chi-1)+1, IM,:,:),p.N_K_chi,mat.Ny,mat.Nz);
            Vsub_mod((Nchi-1)*p.N_K_chi+1:Nchi*p.N_K_chi,:,:) = V_temp;

            m_temp = reshape(m((Nchi-1)*(p.N_K_chi-1)+1:(Nchi)*(p.N_K_chi-1)+1),p.N_K_chi,1);
            m_mod((Nchi-1)*p.N_K_chi+1:Nchi*p.N_K_chi) = m_temp;
        end

        m_mod = m_mod.';

        %%
        
        if (p.do1D == true)
            %% Steady State Case
            if (p.doSteadyState == true)
                [A,rhs] = get_SysM(mat, p, squeeze(Esub(IV, :, IM)));
                

                % if (mat.dg.params.MPC ~= 1)   NOT NECESSARY CAUSE THIS IS
                % DONE IN TRANSIENT FUNCTION
                %     A=mat.dg.params.MPC*A;
                %     rhs=mat.dg.params.MPC*rhs;
                % end
        
        
                rho = A\rhs;
        
                if (mat.dg.params.st_init == true)
                    Atemp{IM} = A;
                    rhstemp{IM} = rhs;
                    rhotemp{IM} = rho;                
                end
                

                %% Calculate carrier density and current density

                %rho_u = (reshape(rho, p.N_chi*p.N_K_chi, p.Nk)*p.phi.').';        %try with only transpose .'
                
                if mat.dg.params.permute_shell == false
                    rho_u = (reshape(rho, p.N_chi*p.N_K_chi, p.Nk)*p.phi.').';     %try with only transpose .'
                elseif mat.dg.params.permute_shell == true
                    rho_u = ((reshape(rho, p.Nk, p.N_chi*p.N_K_chi).')*p.phi.').';
                end
                
                %figure(444),mesh(p.chi_Koord,p.xi_Rechengebiet_Elem,real(rho_u),FaceColor='interp')

                n_DG = (real(rho_u(p.N_xi/2+1,:)+ real(rho_u(p.N_xi/2,:)))/2)';
                j_mode(IV,IM,:) = (p.q*hb./m_mod.*imag(rho_u(p.N_xi/2+1,:)-rho_u(p.N_xi/2,:))/p.delta_xi)';
                %jT      = jT+j_DG_mode;
                
                %Test_n_DG = (squeeze(n_DG(1,1)).*abs(Vsub_mod(1,:,:)).^2/mat.dy/1E-9/mat.dz/1E-9);
                             

                for IN=1:IN_max
                     n_mode(IV,IM,IN,:,:) = (squeeze(n_DG(IN,1)).*squeeze(abs(Vsub_mod(IN,:,:))).^2/mat.dy/1E-9/mat.dz/1E-9);
                end

                RHO_DG(IV,IM,:,:) = rho_u.';
                


                
                %% Plotting
                % if (p.doFV == true)
                %     rho_u = (reshape(rho, p.N_chi*p.N_K_chi, p.N_xi)*p.phi.')';
                %     figure('name', 'Real part of Density Matrix')
                %     mesh(p.chi_Koord*1e9, p.xi_Rechengebiet*1e9,real(rho_u), FaceColor="interp")
                % 
                %     figure('name', 'Imaginary part of Density Matrix')
                %     mesh(p.chi_Koord*1e9, p.xi_Rechengebiet*1e9,imag(rho_u), FaceColor="interp")
                % elseif (p.doDG == true)
                %     rho_u = (reshape(rho, p.N_chi*p.N_K_chi, p.N_xi*p.N_K_xi)*p.phi.')';
                %     figure('name', 'Real part of Density Matrix')
                %     mesh(p.chi_Koord*1e9,p.xi_Koord*1e9,real(rho_u), FaceColor='interp')
                % 
                %     figure('name', 'Imaginary part of Density Matrix')
                %     mesh(p.chi_Koord*1e9,p.xi_Koord*1e9,imag(rho_u), FaceColor='interp')
                % elseif (p.doFEM == true)
                %     rho_u = (reshape(rho, p.N_chi*p.N_K_chi, p.N_xi*(p.N_K_xi-1)+1)*p.phi.')';
                %     figure('name', 'Real part of Density Matrix')
                %     mesh(p.xi_Rechengebiet*1e9,p.chi_Koord*1e9,real(rho_u'), FaceColor='interp')
                % 
                %     figure('name', 'Imaginary part of Density Matrix')
                %     mesh(p.xi_Rechengebiet*1e9,p.chi_Koord*1e9,imag(rho_u'), FaceColor='interp')
                % end
            end
            
           
            
            %% Transient
            % if (p.doTransient == true)
            %     [A,rhs] = get_SysM(p);          % Get Stationary data
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
            %     %% Plotting
            %     if (p.doFV == true)
            %         figure('name', 'Dichte Imaginär'), mesh(p.xi_Rechengebiet*1e9, p.chi_Koord*1e9, imag(rho_T), FaceColor="interp");
            %         figure('name', 'Dichte Real'), mesh(p.xi_Rechengebiet*1e9, p.chi_Koord*1e9, real(rho_T), FaceColor="interp");
            %         figure('name', 'Teilchendichte transient'), mesh(p.chi_Koord*1e9,p.time, n_over_t, FaceColor="interp");
            %         figure('name', 'Stromdichte transient'), mesh(p.chi_Koord*1e9, p.time, j_over_t, FaceColor="interp");
            %     elseif (p.doDG == true)
            %         figure('name', 'Dichte Imaginär'), mesh(p.xi_Koord*1e9, p.chi_Koord*1e9, imag(rho_T), FaceColor="interp");
            %         figure('name', 'Dichte Real'), mesh(p.xi_Koord*1e9, p.chi_Koord*1e9, real(rho_T), FaceColor="interp");
            %         figure('name', 'Teilchendichte transient'), mesh(p.time, p.chi_Koord*1e9, n_over_t, FaceColor="interp");
            %         figure('name', 'Stromdichte transient'), mesh(p.time, p.chi_Koord*1e9, j_over_t, FaceColor="interp");
            %     end
            % end
            
            
        end
    end
end

 nT = squeeze(sum(n_mode, [1 2]));
 jT = squeeze(sum(j_mode, [1 2]));

 if (mat.dg.params.st_init == true)
   DG.A   = Atemp;
   DG.rhs = rhstemp;
   DG.rho = rhotemp;                
 end

 DG.p   = p;
 DG.n   = nT;
 DG.j   = jT;
 DG.RHO = RHO_DG;
 DG.chi = p_top(1,1).chi_Koord;
 DG.xi  = p_top(1,1).xi_Rechengebiet_Elem;





  

