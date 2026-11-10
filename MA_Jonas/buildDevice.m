function mat = buildDevice(device)

% dx = 0.5; % transport direction
% dy = 0.5; % confinement direction

dE = 5E-4;



switch device
    
    case 'DGFET-type-I'
        mat.type = 'DGFET-type-I';
        %n_of_modes = mat.n_of_modes;      

        dx = 0.25; % transport direction
        dy = 0.25; % confinement direction

        n_of_valleys = 1;  % number of valleys
        deg_factor = 1;    % degenarcy factor of the valley
        
        L_c = 10; % channel length
        L_s = 20; % source length
        L_d = 20; % drain length
        
        L_x = L_c + L_s + L_d;
        
        W_c  = 3.0; % channel width
        W_ox = 1.0; % oxid width
        
        L_y = W_c + W_ox + W_ox;
        
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
        
        DEV(1).Eg = Eg_c;
        DEV(1).Eps = Eps_c;
        DEV(1).me_x = me_x_ch;
        DEV(1).me_y = me_y_ch;
        DEV(1).me_z = me_z_ch;
        DEV(1).V = Eg_c;
        DEV(1).Nd = N_s;
        DEV(1).p(1).coord = [0 0];
        DEV(1).p(2).coord = [L_s 0];
        DEV(1).p(3).coord = [L_s W_c];
        DEV(1).p(4).coord = [0 W_c];
        
        DEV(2).Eg = Eg_c;
        DEV(2).Eps = Eps_c;
        DEV(2).me_x = me_x_ch;
        DEV(2).me_y = me_y_ch;
        DEV(2).me_z = me_z_ch;
        DEV(2).V = Eg_c;
        DEV(2).Nd = 0;
        DEV(2).p(1).coord = [L_s 0];
        DEV(2).p(2).coord = [L_s+L_c 0];
        DEV(2).p(3).coord = [L_s+L_c W_c];
        DEV(2).p(4).coord = [L_s W_c];
        
        DEV(3).Eg = Eg_c;
        DEV(3).Eps = Eps_c;
        DEV(3).me_x = me_x_ch;
        DEV(3).me_y = me_y_ch;
        DEV(3).me_z = me_z_ch;
        DEV(3).V = Eg_c;
        DEV(3).Nd = N_d;
        DEV(3).p(1).coord = [L_s+L_c 0];
        DEV(3).p(2).coord = [L_s+L_c+L_d 0];
        DEV(3).p(3).coord = [L_s+L_c+L_d W_c];
        DEV(3).p(4).coord = [L_s+L_c W_c];
        
        %Oxide
        DEV(4).Eg = Eg_ox;
        DEV(4).Eps = Eps_ox;
        DEV(4).me_x = me_x_ox;
        DEV(4).me_y = me_y_ox;
        DEV(4).me_z = me_z_ox;
        DEV(4).V = Eg_c+dEc*(Eg_ox-Eg_c);
        DEV(4).Nd = 0;
        DEV(4).p(1).coord = [0 -W_ox];
        DEV(4).p(2).coord = [L_x -W_ox];
        DEV(4).p(3).coord = [L_x 0];
        DEV(4).p(4).coord = [0 0];
        
        DEV(5).Eg = Eg_ox;
        DEV(5).Eps = Eps_ox;
        DEV(5).me_x = me_x_ox;
        DEV(5).me_y = me_y_ox;
        DEV(5).me_z = me_z_ox;
        DEV(5).V = Eg_c+dEc*(Eg_ox-Eg_c);
        DEV(5).Nd = 0;
        DEV(5).p(1).coord = [0 W_c];
        DEV(5).p(2).coord = [L_x W_c];
        DEV(5).p(3).coord = [L_x W_c+W_ox];
        DEV(5).p(4).coord = [0 W_c+W_c];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MAPPING THE DEVICE TABLE ONTO THE COMPUTATIONAL GRID
        
        mat.x = 0 : dx : L_x;
        mat.y = (0 : dy : L_y)-W_ox;
        
        mat.Nx = length(mat.x);
        mat.Ny = length(mat.y);
        
        mat.V = zeros(mat.Nx, mat.Ny);
        mat.Eg = zeros(mat.Nx, mat.Ny);
        mat.me_x = zeros(n_of_valleys, mat.Nx, mat.Ny);
        mat.me_y = zeros(n_of_valleys, mat.Nx, mat.Ny);
        mat.me_z = zeros(n_of_valleys, mat.Nx, mat.Ny);
        mat.Eps = zeros(mat.Nx, mat.Ny);
        mat.Nd = zeros(mat.Nx, mat.Ny);
        
        
        for IX = 1 : mat.Nx
            for IY = 1 : mat.Ny
                for ID = 1 : length(DEV)
                    
                    if mat.x(IX) >= DEV(ID).p(1).coord(1) && mat.x(IX) <= DEV(ID).p(2).coord(1) && ...
                            mat.y(IY) >=DEV(ID).p(2).coord(2) && mat.y(IY) <= DEV(ID).p(3).coord(2)
                        
                        mat.Eg(IX, IY) = DEV(ID).Eg;
                        mat.V(IX, IY) = DEV(ID).V;
                        mat.Eps(IX, IY) = DEV(ID).Eps;
                        mat.Nd(IX, IY) = DEV(ID).Nd;
                        
                        for IV = 1 : n_of_valleys
                            
                            mat.me_x(IV, IX, IY) = DEV(ID).me_x;
                            mat.me_y(IV, IX, IY) = DEV(ID).me_y;
                            mat.me_z(IV, IX, IY) = DEV(ID).me_z;
                            
                        end     
                        
                    end
                    
                end
            end
        end
        
        mat.L_c = L_c;
        mat.L_s = L_s;
        mat.L_d = L_d;
        
        mat.W_c = W_c;
        mat.W_ox = W_ox;
        
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
        
        mat.dE = dE;
        mat.Temp = 300;
        %mat.n_of_modes = n_of_modes;
        mat.n_of_valleys = n_of_valleys;
        mat.deg_factor = deg_factor;
        
        mat.dx = dx;
        mat.dy = dy;
        
        mat.L_x = L_x;
        
        % boundary assignment for poissons equation
        
        mat.boundary = zeros(mat.Nx, mat.Ny);
        
        gate_idx = mat.x>L_s & mat.x<L_s+L_c;
        
        mat.boundary(gate_idx, 1) = 1;
        mat.boundary(gate_idx, mat.Ny) = 1;
        
    case 'Si'
        
        
    case 'DGFET-GaAs'
        n_of_modes = 1;
        dy = 0.5;
        n_of_valleys = 1;  % number of valleys
        deg_factor = 1;    % degenarcy factor of the valley
        
        L_c = 8; % channel length
        L_s = 15; % source length
        L_d = 15; % drain length
        
        L_x = L_c + L_s + L_d;
        
        W_c  = 3; % channel width
        W_ox = 1; % oxid width
        
        L_y = W_c + W_ox + W_ox;
        
        N_s = 1E+20; % source doping [cm^-3]
        N_d = 1E+20; % drain doping [cm^-3];
        
        Eg_c = 1.424; % bandgap of GaAs
        Eps_c = 12.9;
        
        me_x_ch = 0.063;
        me_y_ch = 0.063;
        me_z_ch = 0.063;
        
        Eg_ox = 8.8;
        Eps_ox = 3.9;
        
        me_x_ox = 0.5;
        me_y_ox = 0.5;
        me_z_ox = 0.5;
        
        Xi = 4.07; % affinity of the channel GaInAs
        phi_m_g = 4.74; % metal work function of the Ag-gate contact
        
        dEc = 0.4;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DEVICE TABLE
        
        DEV(1).Eg = Eg_c;
        DEV(1).Eps = Eps_c;
        DEV(1).me_x = me_x_ch;
        DEV(1).me_y = me_y_ch;
        DEV(1).me_z = me_z_ch;
        DEV(1).V = Eg_c;
        DEV(1).Nd = N_s;
        DEV(1).p(1).coord = [0 0];
        DEV(1).p(2).coord = [L_s 0];
        DEV(1).p(3).coord = [L_s W_c];
        DEV(1).p(4).coord = [0 W_c];
        
        DEV(2).Eg = Eg_c;
        DEV(2).Eps = Eps_c;
        DEV(2).me_x = me_x_ch;
        DEV(2).me_y = me_y_ch;
        DEV(2).me_z = me_z_ch;
        DEV(2).V = Eg_c;
        DEV(2).Nd = 0;
        DEV(2).p(1).coord = [L_s 0];
        DEV(2).p(2).coord = [L_s+L_c 0];
        DEV(2).p(3).coord = [L_s+L_c W_c];
        DEV(2).p(4).coord = [L_s W_c];
        
        DEV(3).Eg = Eg_c;
        DEV(3).Eps = Eps_c;
        DEV(3).me_x = me_x_ch;
        DEV(3).me_y = me_y_ch;
        DEV(3).me_z = me_z_ch;
        DEV(3).V = Eg_c;
        DEV(3).Nd = N_d;
        DEV(3).p(1).coord = [L_s+L_c 0];
        DEV(3).p(2).coord = [L_s+L_c+L_d 0];
        DEV(3).p(3).coord = [L_s+L_c+L_d W_c];
        DEV(3).p(4).coord = [L_s+L_c W_c];
        
        DEV(4).Eg = Eg_ox;
        DEV(4).Eps = Eps_ox;
        DEV(4).me_x = me_x_ox;
        DEV(4).me_y = me_y_ox;
        DEV(4).me_z = me_z_ox;
        DEV(4).V = Eg_c+dEc*(Eg_ox-Eg_c);
        DEV(4).Nd = 0;
        DEV(4).p(1).coord = [0 -W_ox];
        DEV(4).p(2).coord = [L_x -W_ox];
        DEV(4).p(3).coord = [L_x 0];
        DEV(4).p(4).coord = [0 0];
        
        DEV(5).Eg = Eg_ox;
        DEV(5).Eps = Eps_ox;
        DEV(5).me_x = me_x_ox;
        DEV(5).me_y = me_y_ox;
        DEV(5).me_z = me_z_ox;
        DEV(5).V = Eg_c+dEc*(Eg_ox-Eg_c);
        DEV(5).Nd = 0;
        DEV(5).p(1).coord = [0 W_c];
        DEV(5).p(2).coord = [L_x W_c];
        DEV(5).p(3).coord = [L_x W_c+W_ox];
        DEV(5).p(4).coord = [0 W_c+W_c];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MAPPING THE DEVICE TABLE ONTO THE COMPUTATIONAL GRID
        
        mat.x = 0 : dx : L_x;
        mat.y = (0 : dy : L_y)-W_ox;
        
        mat.Nx = length(mat.x);
        mat.Ny = length(mat.y);
        
        mat.V = zeros(mat.Nx, mat.Ny);
        mat.Eg = zeros(mat.Nx, mat.Ny);
        mat.me_x = zeros(n_of_valleys, mat.Nx, mat.Ny);
        mat.me_y = zeros(n_of_valleys, mat.Nx, mat.Ny);
        mat.me_z = zeros(n_of_valleys, mat.Nx, mat.Ny);
        mat.Eps = zeros(mat.Nx, mat.Ny);
        mat.Nd = zeros(mat.Nx, mat.Ny);
        
        
        for IX = 1 : mat.Nx
            for IY = 1 : mat.Ny
                for ID = 1 : length(DEV)
                    
                    if mat.x(IX) >= DEV(ID).p(1).coord(1) && mat.x(IX) <= DEV(ID).p(2).coord(1) && ...
                            mat.y(IY) >=DEV(ID).p(2).coord(2) && mat.y(IY) <= DEV(ID).p(3).coord(2)
                        
                        mat.Eg(IX, IY) = DEV(ID).Eg;
                        mat.V(IX, IY) = DEV(ID).V;
                        mat.Eps(IX, IY) = DEV(ID).Eps;
                        mat.Nd(IX, IY) = DEV(ID).Nd;
                        
                        for IV = 1 : n_of_valleys
                            
                            mat.me_x(IV, IX, IY) = DEV(ID).me_x;
                            mat.me_y(IV, IX, IY) = DEV(ID).me_y;
                            mat.me_z(IV, IX, IY) = DEV(ID).me_z;
                            
                        end     
                        
                    end
                    
                end
            end
        end
        
        mat.L_c = L_c;
        mat.L_s = L_s;
        mat.L_d = L_d;
        
        mat.W_c = W_c;
        mat.W_ox = W_ox;
        
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
        
        mat.dE = dE;
        mat.Temp = 300;
        mat.n_of_modes = n_of_modes;
        mat.n_of_valleys = n_of_valleys;
        mat.deg_factor = deg_factor;
        
        mat.dx = dx;
        mat.dy = dy;
        
        mat.L_x = L_x;
        
        % boundary assignment for poissons equation
        
        mat.boundary = zeros(mat.Nx, mat.Ny);
        
        gate_idx = mat.x>L_s & mat.x<L_s+L_c;
        
        mat.boundary(gate_idx, 1) = 1;
        mat.boundary(gate_idx, mat.Ny) = 1;

    case 'dev1'
        
        n_of_modes = 1;
        n_of_valleys = 1;  % number of valleys
        deg_factor = 1;    % degenarcy factor of the valley
        
        L=zeros(7,1);

        L(1) = 30;
        L(2) = 3;
        L(3) = 3;
        L(4) = 4;
        L(5) = 3;
        L(6) = 3;
        L(7) = 30;


        L_x = sum(L);
        
        W_c  = 5; % channel width
        W_ox = 1; % oxid width
        
        L_y = W_c + W_ox + W_ox;
        
        N_s = 2E+18; % source doping [cm^-3]
        N_d = 2E+18; % drain doping [cm^-3];
        N_ch = 0;1E+15; % channel doping [cm^-3];
        
        Eg_GaAs    = 1.424;
        Eps_GaAs   = 12.9;
        Eg_AlGaAs  = 1.7013;
        Eps_AlGaAs = 12.05;

        me_x_GaAs = 0.063;
        me_y_GaAs = 0.063;
        me_z_GaAs = 0.063;

        me_x_AlGaAs = 0.088;
        me_y_AlGaAs = 0.088;
        me_z_AlGaAs = 0.088;
        
        Eg_ox = 8.8;
        Eps_ox = 3.9;
        
        me_x_ox = 0.5;
        me_y_ox = 0.5;
        me_z_ox = 0.5;

        mat.me_x_ch = 0.063;
        mat.me_y_ch = 0.063;
        mat.me_z_ch = 0.063;
        
        mat.me_x_ox = 0.5;
        mat.me_y_ox = 0.5;
        mat.me_z_ox = 0.5;
        
        Xi = 4.07; % affinity of the channel GaAs
        phi_m_g = 4.74; % metal work function of the Ag-gate contact
        
        dEc = 0.4;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DEVICE TABLE
        
        DEV(1).Eg = Eg_GaAs;
        DEV(1).Eps = Eps_GaAs;
        DEV(1).me_x = me_x_GaAs;
        DEV(1).me_y = me_y_GaAs;
        DEV(1).me_z = me_z_GaAs;
        DEV(1).V = Eg_GaAs;
        DEV(1).Nd = N_s;
        DEV(1).p(1).coord = [0 0];
        DEV(1).p(2).coord = [L(1) 0];
        DEV(1).p(3).coord = [L(1) W_c];
        DEV(1).p(4).coord = [0 W_c];
        
        DEV(2).Eg = Eg_GaAs;
        DEV(2).Eps = Eps_GaAs;
        DEV(2).me_x = me_x_GaAs;
        DEV(2).me_y = me_y_GaAs;
        DEV(2).me_z = me_z_GaAs;
        DEV(2).V = Eg_GaAs;
        DEV(2).Nd = N_ch;
        DEV(2).p(1).coord = [L(1) 0];
        DEV(2).p(2).coord = [L(1)+L(2) 0];
        DEV(2).p(3).coord = [L(1)+L(2) W_c];
        DEV(2).p(4).coord = [L(1) W_c];
        
        DEV(3).Eg = Eg_AlGaAs;
        DEV(3).Eps = Eps_AlGaAs;
        DEV(3).me_x = me_x_AlGaAs;
        DEV(3).me_y = me_y_AlGaAs;
        DEV(3).me_z = me_z_AlGaAs;
        DEV(3).V = Eg_AlGaAs;
        DEV(3).Nd = N_ch;
        DEV(3).p(1).coord = [L(1)+L(2) 0];
        DEV(3).p(2).coord = [L(1)+L(2)+L(3) 0];
        DEV(3).p(3).coord = [L(1)+L(2)+L(3) W_c];
        DEV(3).p(4).coord = [L(1)+L(2) W_c];

        DEV(4).Eg = Eg_GaAs;
        DEV(4).Eps = Eps_GaAs;
        DEV(4).me_x = me_x_GaAs;
        DEV(4).me_y = me_y_GaAs;
        DEV(4).me_z = me_z_GaAs;
        DEV(4).V = Eg_GaAs;
        DEV(4).Nd = N_ch;
        DEV(4).p(1).coord = [L(1)+L(2)+L(3) 0];
        DEV(4).p(2).coord = [L(1)+L(2)+L(3)+L(4) 0];
        DEV(4).p(3).coord = [L(1)+L(2)+L(3)+L(4) W_c];
        DEV(4).p(4).coord = [L(1)+L(2)+L(3) W_c];
        
        DEV(5).Eg = Eg_AlGaAs;
        DEV(5).Eps = Eps_AlGaAs;
        DEV(5).me_x = me_x_AlGaAs;
        DEV(5).me_y = me_y_AlGaAs;
        DEV(5).me_z = me_z_AlGaAs;
        DEV(5).V = Eg_AlGaAs;
        DEV(5).Nd = N_ch;
        DEV(5).p(1).coord = [L(1)+L(2)+L(3)+L(4) 0];
        DEV(5).p(2).coord = [L(1)+L(2)+L(3)+L(4)+L(5) 0];
        DEV(5).p(3).coord = [L(1)+L(2)+L(3)+L(4)+L(5) W_c];
        DEV(5).p(4).coord = [L(1)+L(2)+L(3)+L(4) W_c];

        
        DEV(6).Eg = Eg_GaAs;
        DEV(6).Eps = Eps_GaAs;
        DEV(6).me_x = me_x_GaAs;
        DEV(6).me_y = me_y_GaAs;
        DEV(6).me_z = me_z_GaAs;
        DEV(6).V = Eg_GaAs;
        DEV(6).Nd = N_ch;
        DEV(6).p(1).coord = [L(1)+L(2)+L(3)+L(4)+L(5) 0];
        DEV(6).p(2).coord = [L(1)+L(2)+L(3)+L(4)+L(5)+L(6) 0];
        DEV(6).p(3).coord = [L(1)+L(2)+L(3)+L(4)+L(5)+L(6) W_c];
        DEV(6).p(4).coord = [L(1)+L(2)+L(3)+L(4)+L(5) W_c];

        DEV(7).Eg = Eg_GaAs;
        DEV(7).Eps = Eps_GaAs;
        DEV(7).me_x = me_x_GaAs;
        DEV(7).me_y = me_y_GaAs;
        DEV(7).me_z = me_z_GaAs;
        DEV(7).V = Eg_GaAs;
        DEV(7).Nd = N_d;
        DEV(7).p(1).coord = [L(1)+L(2)+L(3)+L(4)+L(5)+L(6) 0];
        DEV(7).p(2).coord = [L(1)+L(2)+L(3)+L(4)+L(5)+L(6)+L(7) 0];
        DEV(7).p(3).coord = [L(1)+L(2)+L(3)+L(4)+L(5)+L(6)+L(7) W_c];
        DEV(7).p(4).coord = [L(1)+L(2)+L(3)+L(4)+L(5)+L(6) W_c];

        
        DEV(8).Eg = Eg_ox;
        DEV(8).Eps = Eps_ox;
        DEV(8).me_x = me_x_ox;
        DEV(8).me_y = me_y_ox;
        DEV(8).me_z = me_z_ox;
        DEV(8).V = Eg_GaAs+dEc*(Eg_ox-Eg_GaAs);
        DEV(8).Nd = 0;
        DEV(8).p(1).coord = [0 -W_ox];
        DEV(8).p(2).coord = [L_x -W_ox];
        DEV(8).p(3).coord = [L_x 0];
        DEV(8).p(4).coord = [0 0];
        
        DEV(9).Eg = Eg_ox;
        DEV(9).Eps = Eps_ox;
        DEV(9).me_x = me_x_ox;
        DEV(9).me_y = me_y_ox;
        DEV(9).me_z = me_z_ox;
        DEV(9).V = Eg_GaAs+dEc*(Eg_ox-Eg_GaAs);
        DEV(9).Nd = 0;
        DEV(9).p(1).coord = [0 W_c];
        DEV(9).p(2).coord = [L_x W_c];
        DEV(9).p(3).coord = [L_x W_c+W_ox];
        DEV(9).p(4).coord = [0 W_c+W_c];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MAPPING THE DEVICE TABLE ONTO THE COMPUTATIONAL GRID
        
        mat.x = 0 : dx : L_x;
        mat.y = (0 : dy : L_y)-W_ox;
        
        mat.Nx = length(mat.x);
        mat.Ny = length(mat.y);
        
        mat.Id = zeros(mat.Nx, mat.Ny);
        mat.V = zeros(mat.Nx, mat.Ny);
        mat.Eg = zeros(mat.Nx, mat.Ny);
        mat.me_x = zeros(n_of_valleys, mat.Nx, mat.Ny);
        mat.me_y = zeros(n_of_valleys, mat.Nx, mat.Ny);
        mat.me_z = zeros(n_of_valleys, mat.Nx, mat.Ny);
        mat.Eps = zeros(mat.Nx, mat.Ny);
        mat.Nd = zeros(mat.Nx, mat.Ny);
        
        
        for IX = 1 : mat.Nx
            for IY = 1 : mat.Ny
                for ID = 1 : length(DEV)
                    
                    if mat.x(IX) >= DEV(ID).p(1).coord(1) && mat.x(IX) <= DEV(ID).p(2).coord(1) && ...
                            mat.y(IY) >=DEV(ID).p(2).coord(2) && mat.y(IY) <= DEV(ID).p(3).coord(2)
                        
                        mat.Eg(IX, IY) = DEV(ID).Eg;
                        mat.V(IX, IY) = DEV(ID).V;
                        mat.Eps(IX, IY) = DEV(ID).Eps;
                        mat.Nd(IX, IY) = DEV(ID).Nd;
                        mat.Id(IX,IY) = ID;
                        
                        for IV = 1 : n_of_valleys
                            
                            mat.me_x(IV, IX, IY) = DEV(ID).me_x;
                            mat.me_y(IV, IX, IY) = DEV(ID).me_y;
                            mat.me_z(IV, IX, IY) = DEV(ID).me_z;
                            
                        end     
                        
                    end
                    
                end
            end
        end
        
%         mat.L_c = L_c;
%         mat.L_s = L_s;
%         mat.L_d = L_d;
        
        mat.W_c = W_c;
        mat.W_ox = W_ox;
        
%         mat.Eg_c = Eg_c;
%         mat.Eps_c = Eps_c;
        
        mat.Eg_ox = Eg_ox;
        mat.Eps_ox = Eps_ox;
        
        
        mat.Xi = Xi;
        mat.phi_m_g = phi_m_g;
        
        mat.N_s = N_s;
        mat.N_d = N_d;
        
        mat.dE = dE;
        mat.Temp = 77;
        mat.n_of_modes = n_of_modes;
        mat.n_of_valleys = n_of_valleys;
        mat.deg_factor = deg_factor;
        
        mat.dx = dx;
        mat.dy = dy;
        
        mat.L_x = L_x;
        
        % boundary assignment for poissons equation
        
        mat.boundary = zeros(mat.Nx, mat.Ny);
        
        gate_idx = mat.x>(L(1)) & mat.x<(L(1)+L(2)+L(3)+L(4)+L(5)+L(6));
        
%         mat.boundary(gate_idx, 1) = 1;
%         mat.boundary(gate_idx, mat.Ny) = 1;


    case 'dev1g'
        
        n_of_modes = 1;
        n_of_valleys = 1;  % number of valleys
        deg_factor = 1;    % degenarcy factor of the valley
        
        L=zeros(7,1);

        L(1) = 30;
        L(2) = 3;
        L(3) = 3;
        L(4) = 4;
        L(5) = 3;
        L(6) = 3;
        L(7) = 30;


        L_x = sum(L);
        
        W_c  = 5; % channel width
        W_ox = 1; % oxid width
        
        L_y = W_c + W_ox + W_ox;
        
        N_s = 2E+18; % source doping [cm^-3]
        N_d = 2E+18; % drain doping [cm^-3];
        N_ch = 0;1E+15; % channel doping [cm^-3];
        
        Eg_GaAs    = 1.424;
        Eps_GaAs   = 12.9;
        Eg_AlGaAs  = 1.7013;
        Eps_AlGaAs = 12.05;

        me_x_GaAs = 0.063;
        me_y_GaAs = 0.063;
        me_z_GaAs = 0.063;

        me_x_AlGaAs = 0.088;
        me_y_AlGaAs = 0.088;
        me_z_AlGaAs = 0.088;
        
        Eg_ox = 8.8;
        Eps_ox = 3.9;
        
        me_x_ox = 0.5;
        me_y_ox = 0.5;
        me_z_ox = 0.5;

        mat.me_x_ch = 0.063;
        mat.me_y_ch = 0.063;
        mat.me_z_ch = 0.063;
        
        mat.me_x_ox = 0.5;
        mat.me_y_ox = 0.5;
        mat.me_z_ox = 0.5;
        
        Xi = 4.07; % affinity of the channel GaAs
        phi_m_g = 4.74; % metal work function of the Ag-gate contact
        
        dEc = 0.4;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DEVICE TABLE
        
        DEV(1).Eg = Eg_GaAs;
        DEV(1).Eps = Eps_GaAs;
        DEV(1).me_x = me_x_GaAs;
        DEV(1).me_y = me_y_GaAs;
        DEV(1).me_z = me_z_GaAs;
        DEV(1).V = Eg_GaAs;
        DEV(1).Nd = N_s;
        DEV(1).p(1).coord = [0 0];
        DEV(1).p(2).coord = [L(1) 0];
        DEV(1).p(3).coord = [L(1) W_c];
        DEV(1).p(4).coord = [0 W_c];
        
        DEV(2).Eg = Eg_GaAs;
        DEV(2).Eps = Eps_GaAs;
        DEV(2).me_x = me_x_GaAs;
        DEV(2).me_y = me_y_GaAs;
        DEV(2).me_z = me_z_GaAs;
        DEV(2).V = Eg_GaAs;
        DEV(2).Nd = N_ch;
        DEV(2).p(1).coord = [L(1) 0];
        DEV(2).p(2).coord = [L(1)+L(2) 0];
        DEV(2).p(3).coord = [L(1)+L(2) W_c];
        DEV(2).p(4).coord = [L(1) W_c];
        
        DEV(3).Eg = Eg_AlGaAs;
        DEV(3).Eps = Eps_AlGaAs;
        DEV(3).me_x = me_x_AlGaAs;
        DEV(3).me_y = me_y_AlGaAs;
        DEV(3).me_z = me_z_AlGaAs;
        DEV(3).V = Eg_AlGaAs;
        DEV(3).Nd = N_ch;
        DEV(3).p(1).coord = [L(1)+L(2) 0];
        DEV(3).p(2).coord = [L(1)+L(2)+L(3) 0];
        DEV(3).p(3).coord = [L(1)+L(2)+L(3) W_c];
        DEV(3).p(4).coord = [L(1)+L(2) W_c];

        DEV(4).Eg = Eg_GaAs;
        DEV(4).Eps = Eps_GaAs;
        DEV(4).me_x = me_x_GaAs;
        DEV(4).me_y = me_y_GaAs;
        DEV(4).me_z = me_z_GaAs;
        DEV(4).V = Eg_GaAs;
        DEV(4).Nd = N_ch;
        DEV(4).p(1).coord = [L(1)+L(2)+L(3) 0];
        DEV(4).p(2).coord = [L(1)+L(2)+L(3)+L(4) 0];
        DEV(4).p(3).coord = [L(1)+L(2)+L(3)+L(4) W_c];
        DEV(4).p(4).coord = [L(1)+L(2)+L(3) W_c];
        
        DEV(5).Eg = Eg_AlGaAs;
        DEV(5).Eps = Eps_AlGaAs;
        DEV(5).me_x = me_x_AlGaAs;
        DEV(5).me_y = me_y_AlGaAs;
        DEV(5).me_z = me_z_AlGaAs;
        DEV(5).V = Eg_AlGaAs;
        DEV(5).Nd = N_ch;
        DEV(5).p(1).coord = [L(1)+L(2)+L(3)+L(4) 0];
        DEV(5).p(2).coord = [L(1)+L(2)+L(3)+L(4)+L(5) 0];
        DEV(5).p(3).coord = [L(1)+L(2)+L(3)+L(4)+L(5) W_c];
        DEV(5).p(4).coord = [L(1)+L(2)+L(3)+L(4) W_c];

        
        DEV(6).Eg = Eg_GaAs;
        DEV(6).Eps = Eps_GaAs;
        DEV(6).me_x = me_x_GaAs;
        DEV(6).me_y = me_y_GaAs;
        DEV(6).me_z = me_z_GaAs;
        DEV(6).V = Eg_GaAs;
        DEV(6).Nd = N_ch;
        DEV(6).p(1).coord = [L(1)+L(2)+L(3)+L(4)+L(5) 0];
        DEV(6).p(2).coord = [L(1)+L(2)+L(3)+L(4)+L(5)+L(6) 0];
        DEV(6).p(3).coord = [L(1)+L(2)+L(3)+L(4)+L(5)+L(6) W_c];
        DEV(6).p(4).coord = [L(1)+L(2)+L(3)+L(4)+L(5) W_c];

        DEV(7).Eg = Eg_GaAs;
        DEV(7).Eps = Eps_GaAs;
        DEV(7).me_x = me_x_GaAs;
        DEV(7).me_y = me_y_GaAs;
        DEV(7).me_z = me_z_GaAs;
        DEV(7).V = Eg_GaAs;
        DEV(7).Nd = N_d;
        DEV(7).p(1).coord = [L(1)+L(2)+L(3)+L(4)+L(5)+L(6) 0];
        DEV(7).p(2).coord = [L(1)+L(2)+L(3)+L(4)+L(5)+L(6)+L(7) 0];
        DEV(7).p(3).coord = [L(1)+L(2)+L(3)+L(4)+L(5)+L(6)+L(7) W_c];
        DEV(7).p(4).coord = [L(1)+L(2)+L(3)+L(4)+L(5)+L(6) W_c];

        
        DEV(8).Eg = Eg_ox;
        DEV(8).Eps = Eps_ox;
        DEV(8).me_x = me_x_ox;
        DEV(8).me_y = me_y_ox;
        DEV(8).me_z = me_z_ox;
        DEV(8).V = Eg_GaAs+dEc*(Eg_ox-Eg_GaAs);
        DEV(8).Nd = 0;
        DEV(8).p(1).coord = [0 -W_ox];
        DEV(8).p(2).coord = [L_x -W_ox];
        DEV(8).p(3).coord = [L_x 0];
        DEV(8).p(4).coord = [0 0];
        
        DEV(9).Eg = Eg_ox;
        DEV(9).Eps = Eps_ox;
        DEV(9).me_x = me_x_ox;
        DEV(9).me_y = me_y_ox;
        DEV(9).me_z = me_z_ox;
        DEV(9).V = Eg_GaAs+dEc*(Eg_ox-Eg_GaAs);
        DEV(9).Nd = 0;
        DEV(9).p(1).coord = [0 W_c];
        DEV(9).p(2).coord = [L_x W_c];
        DEV(9).p(3).coord = [L_x W_c+W_ox];
        DEV(9).p(4).coord = [0 W_c+W_c];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MAPPING THE DEVICE TABLE ONTO THE COMPUTATIONAL GRID
        
        mat.x = 0 : dx : L_x;
        mat.y = (0 : dy : L_y)-W_ox;
        
        mat.Nx = length(mat.x);
        mat.Ny = length(mat.y);
        
        mat.Id = zeros(mat.Nx, mat.Ny);
        mat.V = zeros(mat.Nx, mat.Ny);
        mat.Eg = zeros(mat.Nx, mat.Ny);
        mat.me_x = zeros(n_of_valleys, mat.Nx, mat.Ny);
        mat.me_y = zeros(n_of_valleys, mat.Nx, mat.Ny);
        mat.me_z = zeros(n_of_valleys, mat.Nx, mat.Ny);
        mat.Eps = zeros(mat.Nx, mat.Ny);
        mat.Nd = zeros(mat.Nx, mat.Ny);
        
        
        for IX = 1 : mat.Nx
            for IY = 1 : mat.Ny
                for ID = 1 : length(DEV)
                    
                    if mat.x(IX) >= DEV(ID).p(1).coord(1) && mat.x(IX) <= DEV(ID).p(2).coord(1) && ...
                            mat.y(IY) >=DEV(ID).p(2).coord(2) && mat.y(IY) <= DEV(ID).p(3).coord(2)
                        
                        mat.Eg(IX, IY) = DEV(ID).Eg;
                        mat.V(IX, IY) = DEV(ID).V;
                        mat.Eps(IX, IY) = DEV(ID).Eps;
                        mat.Nd(IX, IY) = DEV(ID).Nd;
                        mat.Id(IX,IY) = ID;
                        
                        for IV = 1 : n_of_valleys
                            
                            mat.me_x(IV, IX, IY) = DEV(ID).me_x;
                            mat.me_y(IV, IX, IY) = DEV(ID).me_y;
                            mat.me_z(IV, IX, IY) = DEV(ID).me_z;
                            
                        end     
                        
                    end
                    
                end
            end
        end
        
%         mat.L_c = L_c;
%         mat.L_s = L_s;
%         mat.L_d = L_d;
        
        mat.W_c = W_c;
        mat.W_ox = W_ox;
        
%         mat.Eg_c = Eg_c;
%         mat.Eps_c = Eps_c;
        
        mat.Eg_ox = Eg_ox;
        mat.Eps_ox = Eps_ox;
        
        
        mat.Xi = Xi;
        mat.phi_m_g = phi_m_g;
        
        mat.N_s = N_s;
        mat.N_d = N_d;
        
        mat.dE = dE;
        mat.Temp = 77;
        mat.n_of_modes = n_of_modes;
        mat.n_of_valleys = n_of_valleys;
        mat.deg_factor = deg_factor;
        
        mat.dx = dx;
        mat.dy = dy;
        
        mat.L_x = L_x;
        
        % boundary assignment for poissons equation
        
        mat.boundary = zeros(mat.Nx, mat.Ny);
        
        gate_idx = mat.x>(L(1)) & mat.x<(L(1)+L(2)+L(3)+L(4)+L(5)+L(6));
        
        mat.boundary(gate_idx, 1) = 1;
        mat.boundary(gate_idx, mat.Ny) = 1;




        case 'dev2'
        
        n_of_modes = 3;
        n_of_valleys = 1;  % number of valleys
        deg_factor = 1;    % degenarcy factor of the valley
        
        L=zeros(7,1);

        L(1) = 30;
        L(2) = 3;
        L(3) = 3;
        L(4) = 4;
        L(5) = 3;
        L(6) = 3;
        L(7) = 30;


        L_x = sum(L);
        
        W_c  = 20; % channel width
        W_ox = 1; % oxid width
        
        L_y = W_c + W_ox + W_ox;
        
        N_s = 2E+18; % source doping [cm^-3]
        N_d = 2E+18; % drain doping [cm^-3];
        N_ch = 0;1E+15; % channel doping [cm^-3];
        
        Eg_GaAs    = 1.424;
        Eps_GaAs   = 12.9;
        Eg_AlGaAs  = 1.7013;
        Eps_AlGaAs = 12.05;

        me_x_GaAs = 0.063;
        me_y_GaAs = 0.063;
        me_z_GaAs = 0.063;

        me_x_AlGaAs = 0.088;
        me_y_AlGaAs = 0.088;
        me_z_AlGaAs = 0.088;
        
        Eg_ox = 8.8;
        Eps_ox = 3.9;
        
        me_x_ox = 0.5;
        me_y_ox = 0.5;
        me_z_ox = 0.5;

        mat.me_x_ch = 0.063;
        mat.me_y_ch = 0.063;
        mat.me_z_ch = 0.063;
        
        mat.me_x_ox = 0.5;
        mat.me_y_ox = 0.5;
        mat.me_z_ox = 0.5;
        
        Xi = 4.07; % affinity of the channel GaAs
        phi_m_g = 4.74; % metal work function of the Ag-gate contact
        
        dEc = 0.4;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DEVICE TABLE
        
        DEV(1).Eg = Eg_GaAs;
        DEV(1).Eps = Eps_GaAs;
        DEV(1).me_x = me_x_GaAs;
        DEV(1).me_y = me_y_GaAs;
        DEV(1).me_z = me_z_GaAs;
        DEV(1).V = Eg_GaAs;
        DEV(1).Nd = N_s;
        DEV(1).p(1).coord = [0 0];
        DEV(1).p(2).coord = [L(1) 0];
        DEV(1).p(3).coord = [L(1) W_c];
        DEV(1).p(4).coord = [0 W_c];
        
        DEV(2).Eg = Eg_GaAs;
        DEV(2).Eps = Eps_GaAs;
        DEV(2).me_x = me_x_GaAs;
        DEV(2).me_y = me_y_GaAs;
        DEV(2).me_z = me_z_GaAs;
        DEV(2).V = Eg_GaAs;
        DEV(2).Nd = N_ch;
        DEV(2).p(1).coord = [L(1) 0];
        DEV(2).p(2).coord = [L(1)+L(2) 0];
        DEV(2).p(3).coord = [L(1)+L(2) W_c];
        DEV(2).p(4).coord = [L(1) W_c];
        
        DEV(3).Eg = Eg_AlGaAs;
        DEV(3).Eps = Eps_AlGaAs;
        DEV(3).me_x = me_x_AlGaAs;
        DEV(3).me_y = me_y_AlGaAs;
        DEV(3).me_z = me_z_AlGaAs;
        DEV(3).V = Eg_AlGaAs;
        DEV(3).Nd = N_ch;
        DEV(3).p(1).coord = [L(1)+L(2) 0];
        DEV(3).p(2).coord = [L(1)+L(2)+L(3) 0];
        DEV(3).p(3).coord = [L(1)+L(2)+L(3) W_c];
        DEV(3).p(4).coord = [L(1)+L(2) W_c];

        DEV(4).Eg = Eg_GaAs;
        DEV(4).Eps = Eps_GaAs;
        DEV(4).me_x = me_x_GaAs;
        DEV(4).me_y = me_y_GaAs;
        DEV(4).me_z = me_z_GaAs;
        DEV(4).V = Eg_GaAs;
        DEV(4).Nd = N_ch;
        DEV(4).p(1).coord = [L(1)+L(2)+L(3) 0];
        DEV(4).p(2).coord = [L(1)+L(2)+L(3)+L(4) 0];
        DEV(4).p(3).coord = [L(1)+L(2)+L(3)+L(4) W_c];
        DEV(4).p(4).coord = [L(1)+L(2)+L(3) W_c];
        
        DEV(5).Eg = Eg_AlGaAs;
        DEV(5).Eps = Eps_AlGaAs;
        DEV(5).me_x = me_x_AlGaAs;
        DEV(5).me_y = me_y_AlGaAs;
        DEV(5).me_z = me_z_AlGaAs;
        DEV(5).V = Eg_AlGaAs;
        DEV(5).Nd = N_ch;
        DEV(5).p(1).coord = [L(1)+L(2)+L(3)+L(4) 0];
        DEV(5).p(2).coord = [L(1)+L(2)+L(3)+L(4)+L(5) 0];
        DEV(5).p(3).coord = [L(1)+L(2)+L(3)+L(4)+L(5) W_c];
        DEV(5).p(4).coord = [L(1)+L(2)+L(3)+L(4) W_c];

        
        DEV(6).Eg = Eg_GaAs;
        DEV(6).Eps = Eps_GaAs;
        DEV(6).me_x = me_x_GaAs;
        DEV(6).me_y = me_y_GaAs;
        DEV(6).me_z = me_z_GaAs;
        DEV(6).V = Eg_GaAs;
        DEV(6).Nd = N_ch;
        DEV(6).p(1).coord = [L(1)+L(2)+L(3)+L(4)+L(5) 0];
        DEV(6).p(2).coord = [L(1)+L(2)+L(3)+L(4)+L(5)+L(6) 0];
        DEV(6).p(3).coord = [L(1)+L(2)+L(3)+L(4)+L(5)+L(6) W_c];
        DEV(6).p(4).coord = [L(1)+L(2)+L(3)+L(4)+L(5) W_c];

        DEV(7).Eg = Eg_GaAs;
        DEV(7).Eps = Eps_GaAs;
        DEV(7).me_x = me_x_GaAs;
        DEV(7).me_y = me_y_GaAs;
        DEV(7).me_z = me_z_GaAs;
        DEV(7).V = Eg_GaAs;
        DEV(7).Nd = N_d;
        DEV(7).p(1).coord = [L(1)+L(2)+L(3)+L(4)+L(5)+L(6) 0];
        DEV(7).p(2).coord = [L(1)+L(2)+L(3)+L(4)+L(5)+L(6)+L(7) 0];
        DEV(7).p(3).coord = [L(1)+L(2)+L(3)+L(4)+L(5)+L(6)+L(7) W_c];
        DEV(7).p(4).coord = [L(1)+L(2)+L(3)+L(4)+L(5)+L(6) W_c];

        
        DEV(8).Eg = Eg_ox;
        DEV(8).Eps = Eps_ox;
        DEV(8).me_x = me_x_ox;
        DEV(8).me_y = me_y_ox;
        DEV(8).me_z = me_z_ox;
        DEV(8).V = Eg_GaAs+dEc*(Eg_ox-Eg_GaAs);
        DEV(8).Nd = 0;
        DEV(8).p(1).coord = [0 -W_ox];
        DEV(8).p(2).coord = [L_x -W_ox];
        DEV(8).p(3).coord = [L_x 0];
        DEV(8).p(4).coord = [0 0];
        
        DEV(9).Eg = Eg_ox;
        DEV(9).Eps = Eps_ox;
        DEV(9).me_x = me_x_ox;
        DEV(9).me_y = me_y_ox;
        DEV(9).me_z = me_z_ox;
        DEV(9).V = Eg_GaAs+dEc*(Eg_ox-Eg_GaAs);
        DEV(9).Nd = 0;
        DEV(9).p(1).coord = [0 W_c];
        DEV(9).p(2).coord = [L_x W_c];
        DEV(9).p(3).coord = [L_x W_c+W_ox];
        DEV(9).p(4).coord = [0 W_c+W_c];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MAPPING THE DEVICE TABLE ONTO THE COMPUTATIONAL GRID
        
        mat.x = 0 : dx : L_x;
        mat.y = (0 : dy : L_y)-W_ox;
        
        mat.Nx = length(mat.x);
        mat.Ny = length(mat.y);
        
        mat.Id = zeros(mat.Nx, mat.Ny);
        mat.V = zeros(mat.Nx, mat.Ny);
        mat.Eg = zeros(mat.Nx, mat.Ny);
        mat.me_x = zeros(n_of_valleys, mat.Nx, mat.Ny);
        mat.me_y = zeros(n_of_valleys, mat.Nx, mat.Ny);
        mat.me_z = zeros(n_of_valleys, mat.Nx, mat.Ny);
        mat.Eps = zeros(mat.Nx, mat.Ny);
        mat.Nd = zeros(mat.Nx, mat.Ny);
        
        
        for IX = 1 : mat.Nx
            for IY = 1 : mat.Ny
                for ID = 1 : length(DEV)
                    
                    if mat.x(IX) >= DEV(ID).p(1).coord(1) && mat.x(IX) <= DEV(ID).p(2).coord(1) && ...
                            mat.y(IY) >=DEV(ID).p(2).coord(2) && mat.y(IY) <= DEV(ID).p(3).coord(2)
                        
                        mat.Eg(IX, IY) = DEV(ID).Eg;
                        mat.V(IX, IY) = DEV(ID).V;
                        mat.Eps(IX, IY) = DEV(ID).Eps;
                        mat.Nd(IX, IY) = DEV(ID).Nd;
                        mat.Id(IX,IY) = ID;
                        
                        for IV = 1 : n_of_valleys
                            
                            mat.me_x(IV, IX, IY) = DEV(ID).me_x;
                            mat.me_y(IV, IX, IY) = DEV(ID).me_y;
                            mat.me_z(IV, IX, IY) = DEV(ID).me_z;
                            
                        end     
                        
                    end
                    
                end
            end
        end
        
%         mat.L_c = L_c;
%         mat.L_s = L_s;
%         mat.L_d = L_d;
        
        mat.W_c = W_c;
        mat.W_ox = W_ox;
        
%         mat.Eg_c = Eg_c;
%         mat.Eps_c = Eps_c;
        
        mat.Eg_ox = Eg_ox;
        mat.Eps_ox = Eps_ox;
        
        
        mat.Xi = Xi;
        mat.phi_m_g = phi_m_g;
        
        mat.N_s = N_s;
        mat.N_d = N_d;
        
        mat.dE = dE;
        mat.Temp = 77;
        mat.n_of_modes = n_of_modes;
        mat.n_of_valleys = n_of_valleys;
        mat.deg_factor = deg_factor;
        
        mat.dx = dx;
        mat.dy = dy;
        
        mat.L_x = L_x;
        
        % boundary assignment for poissons equation
        
        mat.boundary = zeros(mat.Nx, mat.Ny);
        
        gate_idx = mat.x>(L(1)+L(2)) & mat.x<(L(1)+L(2)+L(3)+L(4)+L(5));
        
%         mat.boundary(gate_idx, 1) = 1;
%         mat.boundary(gate_idx, mat.Ny) = 1;
end


