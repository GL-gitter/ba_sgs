function mat = buildDeviceCMS_3D(device, arg_dx,arg_dy,arg_dz, arg_w_chy,arg_w_oxy,arg_w_chz,arg_w_oxz)

% dx = 0.5; % transport direction
% dy = 0.5; % confinement direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% device = 'GAAFET'
% dx      = 0.25; %0,25
% dy      = 0.25; %0,25
% w_chy    = 4;
% w_chz = w_chy;
% w_oxy    = 1;
% w_oxz = w_oxy;
%
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dE=1e3;
switch device

    case 'GAAFET'
        mat.type = 'GAAFET';

        dx = arg_dx; % transport direction
        dy = arg_dy; % confinement direction
        dz = arg_dz; % confinement direction

        n_of_valleys = 1;  % number of valleys
        deg_factor = 1;    % degenarcy factor of the valley

        L_c = 10; % channel length
        L_s = 21; % source length
        L_d = 21; % drain length

        L_x = L_c + L_s + L_d;

        W_cy  = arg_w_chy; % channel width y
        W_cz  = arg_w_chz; % channel width z
        W_oxy = arg_w_oxy; % oxid width y
        W_oxz = arg_w_oxz; % oxid width z

        L_y = W_cy + W_oxy + W_oxy;
        L_z = W_cz + W_oxz + W_oxz;

        N_s = 2E+19; % source doping [cm^-3]
        N_d = 2E+19; % drain doping [cm^-3];

        Eg_c = 0.74; % bandgap of InGaAs
        Eps_c = 13.9;

        me_x_ch = 0.041;
        me_y_ch = 0.041;
        me_z_ch = 0.041;

        Eg_ox = 8.8;
        Eps_ox = 3.9;

        me_x_ox = 0.5;
        me_y_ox = 0.5;
        me_z_ox = 0.5;

        Xi = 4.5; % affinity of the channel GaInAs
        phi_m_g = 4.74; % metal work function of the Ag-gate contact

        dEc = 0.4;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DEVICE TABLE

        % OXIDE
        DEV(1).Eg = Eg_ox;
        DEV(1).Eps = Eps_ox;
        DEV(1).me_x = me_x_ox;
        DEV(1).me_y = me_y_ox;
        DEV(1).me_z = me_z_ox;
        DEV(1).V = Eg_c+dEc*(Eg_ox-Eg_c);
        DEV(1).Nd = 0;

        % Channel
        DEV(2).Eg = Eg_c;
        DEV(2).Eps = Eps_c;
        DEV(2).me_x = me_x_ch;
        DEV(2).me_y = me_y_ch;
        DEV(2).me_z = me_z_ch;
        DEV(2).V = Eg_c;
        DEV(2).Nd = 0;                     %% Changed for Mathias Paper, normally channel is doped like Source
        DEV(2).p(1).coord = [0 0];
        DEV(2).p(2).coord = [W_cy 0];
        DEV(2).p(3).coord = [W_cy W_cz];
        DEV(2).p(4).coord = [0 W_cz];

        % source / drain
        DEV(3).Eg = Eg_c;
        DEV(3).Eps = Eps_c;
        DEV(3).me_x = me_x_ch;
        DEV(3).me_y = me_y_ch;
        DEV(3).me_z = me_z_ch;
        DEV(3).V = Eg_c;
        DEV(3).Nd = N_s;
        DEV(3).p(1).coord = [0 0];
        DEV(3).p(2).coord = [W_cy 0];
        DEV(3).p(3).coord = [W_cy W_cz];
        DEV(3).p(4).coord = [0 W_cz];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MAPPING THE DEVICE TABLE ONTO THE COMPUTATIONAL GRID

        mat.x = 0 : dx : L_x;
        mat.y = (0 : dy : L_y)-W_oxy;
        mat.z = (0 : dz : L_z)-W_oxz;

        mat.Nx = length(mat.x);
        mat.Ny = length(mat.y);
        mat.Nz = length(mat.z);

        mat.V = zeros(mat.Nx, mat.Ny, mat.Nz);
        mat.Eg = zeros(mat.Nx, mat.Ny, mat.Nz);
        mat.me_x = zeros(n_of_valleys, mat.Nx, mat.Ny, mat.Nz);
        mat.me_y = zeros(n_of_valleys, mat.Nx, mat.Ny, mat.Nz);
        mat.me_z = zeros(n_of_valleys, mat.Nx, mat.Ny, mat.Nz);
        mat.Eps = zeros(mat.Nx, mat.Ny, mat.Nz);
        mat.Nd = zeros(mat.Nx, mat.Ny, mat.Nz);


        %% Changed this snippet for Mathias Paper
        mat.Eg(:,:,:) = DEV(1).Eg;
        mat.V(:,:,:) = DEV(1).V;
        mat.Eps(:,:,:) = DEV(1).Eps;


        mat.me_x(:,:,:) = DEV(1).me_x;
        mat.me_y(:,:,:) = DEV(1).me_y;
        mat.me_z(:,:,:) = DEV(1).me_z;


        % for IX = 1 : mat.Nx
        %     for IY = 1 : mat.Ny
        %         for IZ = 1 : mat.Nz
        % 
        %             if (mat.x(IX)<L_s || mat.x(IX)>L_x-L_d)
        % 
        %                 AREA=3;
        %             else
        %                 AREA=2;
        %             end
        % 
        % 
        %             if mat.y(IY) > DEV(AREA).p(1).coord(1) && mat.y(IY) < DEV(AREA).p(2).coord(1) && ...
        %                     mat.z(IZ) >DEV(AREA).p(2).coord(2) && mat.z(IZ) < DEV(AREA).p(3).coord(2)
        % 
        % 
        %                 %channel
        %                 mat.Eg(IX, IY, IZ) = DEV(AREA).Eg;
        %                 mat.V(IX, IY, IZ) = DEV(AREA).V;
        %                 mat.Eps(IX, IY, IZ) = DEV(AREA).Eps;
        % 
        %                 mat.Nd(IX, IY, IZ) = DEV(AREA).Nd;
        % 
        % 
        %                 for IV = 1 : n_of_valleys
        % 
        %                     mat.me_x(IV, IX, IY, IZ) = DEV(AREA).me_x;
        %                     mat.me_y(IV, IX, IY, IZ) = DEV(AREA).me_y;
        %                     mat.me_z(IV, IX, IY, IZ) = DEV(AREA).me_z;
        % 
        %                 end
        % 
        % 
        %             end
        % 
        %         end
        %     end
        % end

        %%

        %Uncomment this for normal version
        for IX = 1 : mat.Nx
            for IY = 1 : mat.Ny
                for IZ = 1 : mat.Nz
                    if mat.y(IY) > DEV(2).p(1).coord(1) && mat.y(IY) < DEV(2).p(2).coord(1) && ...
                            mat.z(IZ) >DEV(2).p(2).coord(2) && mat.z(IZ) < DEV(2).p(3).coord(2)

                        %channel
                        mat.Eg(IX, IY, IZ) = DEV(2).Eg;
                        mat.V(IX, IY, IZ) = DEV(2).V;
                        mat.Eps(IX, IY, IZ) = DEV(2).Eps;

                        if mat.x(IX)<=L_s
                            mat.Nd(IX, IY, IZ) = N_s;
                        elseif mat.x(IX)>=L_s+L_c
                            mat.Nd(IX, IY, IZ) = N_d;
                        else
                            mat.Nd(IX, IY, IZ) = 0;
                        end

                        for IV = 1 : n_of_valleys

                            mat.me_x(IV, IX, IY, IZ) = DEV(2).me_x;
                            mat.me_y(IV, IX, IY, IZ) = DEV(2).me_y;
                            mat.me_z(IV, IX, IY, IZ) = DEV(2).me_z;

                        end

                    else
                        %oxide

                        mat.Eg(IX, IY, IZ) = DEV(1).Eg;
                        mat.V(IX, IY, IZ) = DEV(1).V;
                        mat.Eps(IX, IY, IZ) = DEV(1).Eps;
                        mat.Nd(IX, IY, IZ) = DEV(1).Nd;

                        for IV = 1 : n_of_valleys
                            mat.me_x(IV, IX, IY, IZ) = DEV(1).me_x;
                            mat.me_y(IV, IX, IY, IZ) = DEV(1).me_y;
                            mat.me_z(IV, IX, IY, IZ) = DEV(1).me_z;
                        end

                    end

                end
            end
        end

        mat.L_c = L_c;
        mat.L_s = L_s;
        mat.L_d = L_d;

        mat.W_cy = W_cy;
        mat.W_cz = W_cz;
        mat.W_oxy = W_oxy;
        mat.W_oxz = W_oxz;

        mat.Eg_c = Eg_c;
        mat.Eps_c = Eps_c;

        mat.Eg_ox = Eg_ox;
        mat.Eps_ox = Eps_ox;

        mat.me_x_ch = me_x_ch;
        mat.me_y_ch = me_y_ch;
        mat.me_z_ch = me_z_ch;

        mat.me_x_ox = me_x_ox;
        mat.me_y_ox = me_y_ox;
        mat.me_z_ox = me_z_ox;

        mat.Xi = Xi;
        mat.phi_m_g = phi_m_g;

        mat.N_s = N_s;
        mat.N_d = N_d;

        %mat.dE = dE;
        mat.Temp = 300;
        %mat.n_of_modes = n_of_modes;
        mat.n_of_valleys = n_of_valleys;
        mat.deg_factor = deg_factor;

        mat.dx = dx;
        mat.dy = dy;
        mat.dz = dz;

        mat.L_x = L_x;

        % boundary assignment for poissons equation

        mat.boundary = zeros(mat.Nx, mat.Ny, mat.Nz);



        %mat.gate_idx = mat.x>L_s & mat.x<L_s+L_c;
        mat.gate_idx = mat.x>L_s & mat.x<L_s+L_c;


        mat.boundary(mat.gate_idx, 1, :) = 1;
        mat.boundary(mat.gate_idx, end, :) = 1;
        mat.boundary(mat.gate_idx, :,1) = 1;
        mat.boundary(mat.gate_idx, :,end) = 1;


end


