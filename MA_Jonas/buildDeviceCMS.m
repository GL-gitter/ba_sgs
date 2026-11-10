function mat = buildDeviceCMS(device, arg_dx,arg_dy,arg_w_ch,arg_w_ox,arg_offsetL, arg_offsetU,arg_taper)

% dx = 0.5; % transport direction
% dy = 0.5; % confinement direction


dE=1e-3;
switch device
    
    case 'DGFET'
        mat.type = 'DGFET';
        dx = arg_dx; % transport direction
        dy = arg_dy; % confinement direction

        n_of_valleys = 1;  % number of valleys
        deg_factor = 1;    % degenarcy factor of the valley
        
        L_c = 10; % channel length
        L_s = 20; % source length
        L_d = 20; % drain length
        
        L_x = L_c + L_s + L_d;
        
        W_c  = arg_w_ch; % channel width
        W_ox = 1; % oxid width
        
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
        % mat.n_of_modes = n_of_modes;
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
        

    case 'SQDGFET'
        % n_of_modes = 1;
        mat.type = 'SQDGFET';
        dx = arg_dx; % transport direction
        dy = arg_dy; % confinement direction

        n_of_valleys = 1;  % number of valleys
        deg_factor = 1;    % degenarcy factor of the valley
        
        L_c = 10; % channel length
        L_s = 20; % source length
        L_d = 20; % drain length
        
        L_x = L_c + L_s + L_d;
        
        W_c  = arg_w_ch; % channel width
        SQU  = arg_offsetU;
        SQL  = arg_offsetL;
        W_ox = 1; % oxid width
        
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
        DEV(2).p(1).coord = [L_s 0+SQL];
        DEV(2).p(2).coord = [L_s+L_c 0+SQL];
        DEV(2).p(3).coord = [L_s+L_c W_c-SQU];
        DEV(2).p(4).coord = [L_s W_c-SQU];
        
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
        DEV(4).p(3).coord = [L_x W_c+W_ox];
        DEV(4).p(4).coord = [0 W_c+W_ox];
        
        % DEV(5).Eg = Eg_ox;
        % DEV(5).Eps = Eps_ox;
        % DEV(5).me_x = me_x_ox;
        % DEV(5).me_y = me_y_ox;
        % DEV(5).me_z = me_z_ox;
        % DEV(5).V = Eg_c+dEc*(Eg_ox-Eg_c);
        % DEV(5).Nd = 0;
        % DEV(5).p(1).coord = [0 W_c];
        % DEV(5).p(2).coord = [L_x W_c];
        % DEV(5).p(3).coord = [L_x W_c+W_ox];
        % DEV(5).p(4).coord = [0 W_c+W_c];
        % 
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
        
        DEVo = [4,3,1,2];
        for IX = 1 : mat.Nx
            for IY = 1 : mat.Ny
                for IDx = 1 : length(DEV)
                    ID = DEVo(IDx);
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
        % mat.n_of_modes = n_of_modes;
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

    case 'TAPSQDGFET'
        % n_of_modes = 1;
        mat.type = 'TAPSQDGFET';
        dx = arg_dx; % transport direction
        dy = arg_dy; % confinement direction

        n_of_valleys = 1;  % number of valleys
        deg_factor = 1;    % degenarcy factor of the valley
        
        L_c = 10; % channel length
        L_s = 20; % source length
        L_d = 20; % drain length
        
        L_x = L_c + L_s + L_d;
        
        W_c  = arg_w_ch; % channel width

        SQU  = arg_offsetU; % squeeze upper
        SQL  = arg_offsetL; % squeeze lower

        taper = arg_taper;

        W_ox = arg_w_ox; % oxid width
        
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
        % me_x_ox = 0.041;
        % me_y_ox = 0.041;
        % me_z_ox = 0.041;
        
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
        DEV(1).p(2).coord = [L_s-taper/2 0];
        DEV(1).p(3).coord = [L_s-taper/2 W_c];
        DEV(1).p(4).coord = [0 W_c];
        
        DEV(2).Eg = Eg_c;
        DEV(2).Eps = Eps_c;
        DEV(2).me_x = me_x_ch;
        DEV(2).me_y = me_y_ch;
        DEV(2).me_z = me_z_ch;
        DEV(2).V = Eg_c;
        DEV(2).Nd = 0;
        DEV(2).p(1).coord = [L_s-taper/2 0+SQL*0];
        DEV(2).p(2).coord = [L_s+taper/2 0+SQL*1];
        DEV(2).p(3).coord = [L_s+taper/2 W_c-SQU*1];
        DEV(2).p(4).coord = [L_s-taper/2 W_c-SQU*0];



        DEV(3).Eg = Eg_c;
        DEV(3).Eps = Eps_c;
        DEV(3).me_x = me_x_ch;
        DEV(3).me_y = me_y_ch;
        DEV(3).me_z = me_z_ch;
        DEV(3).V = Eg_c;
        DEV(3).Nd = 0;
        DEV(3).p(1).coord = [L_s+taper/2 0+SQL];
        DEV(3).p(2).coord = [L_s+L_c-taper/2 0+SQL];
        DEV(3).p(3).coord = [L_s+L_c-taper/2 W_c-SQU];
        DEV(3).p(4).coord = [L_s+taper/2 W_c-SQU];

        DEV(4).Eg = Eg_c;
        DEV(4).Eps = Eps_c;
        DEV(4).me_x = me_x_ch;
        DEV(4).me_y = me_y_ch;
        DEV(4).me_z = me_z_ch;
        DEV(4).V = Eg_c;
        DEV(4).Nd = 0;
        DEV(4).p(1).coord = [L_s+L_c-taper/2 0+SQL*1];
        DEV(4).p(2).coord = [L_s+L_c+taper/2 0+SQL*0];
        DEV(4).p(3).coord = [L_s+L_c+taper/2 W_c-SQU*0];
        DEV(4).p(4).coord = [L_s+L_c-taper/2 W_c-SQU*1];
        
        DEV(5).Eg = Eg_c;
        DEV(5).Eps = Eps_c;
        DEV(5).me_x = me_x_ch;
        DEV(5).me_y = me_y_ch;
        DEV(5).me_z = me_z_ch;
        DEV(5).V = Eg_c;
        DEV(5).Nd = N_d;
        DEV(5).p(1).coord = [L_s+L_c+taper/2 0];
        DEV(5).p(2).coord = [L_s+L_c+L_d 0];
        DEV(5).p(3).coord = [L_s+L_c+L_d W_c];
        DEV(5).p(4).coord = [L_s+L_c+taper/2 W_c];
        
        DEV(6).Eg = Eg_ox;
        DEV(6).Eps = Eps_ox;
        DEV(6).me_x = me_x_ox;
        DEV(6).me_y = me_y_ox;
        DEV(6).me_z = me_z_ox;
        DEV(6).V = Eg_c+dEc*(Eg_ox-Eg_c);
        DEV(6).Nd = 0;
        DEV(6).p(1).coord = [0 -W_ox];
        DEV(6).p(2).coord = [L_x -W_ox];
        DEV(6).p(3).coord = [L_x W_c+W_ox];
        DEV(6).p(4).coord = [0 W_c+W_ox];
        
        % DEV(5).Eg = Eg_ox;
        % DEV(5).Eps = Eps_ox;
        % DEV(5).me_x = me_x_ox;
        % DEV(5).me_y = me_y_ox;
        % DEV(5).me_z = me_z_ox;
        % DEV(5).V = Eg_c+dEc*(Eg_ox-Eg_c);
        % DEV(5).Nd = 0;
        % DEV(5).p(1).coord = [0 W_c];
        % DEV(5).p(2).coord = [L_x W_c];
        % DEV(5).p(3).coord = [L_x W_c+W_ox];
        % DEV(5).p(4).coord = [0 W_c+W_c];
        % 
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
        
        DEVo = [6, 2, 4,3,1,5];
        for IX = 1 : mat.Nx
            for IY = 1 : mat.Ny
                for IDx = 1 : length(DEV)
                    ID = DEVo(IDx);
                    % if mat.x(IX) >= DEV(ID).p(1).coord(1) && mat.x(IX) <= DEV(ID).p(2).coord(1) && ...
                    %         mat.y(IY) >=DEV(ID).p(2).coord(2) && mat.y(IY) <= DEV(ID).p(3).coord(2)
                    if inpolygon(mat.x(IX), mat.y(IY),[DEV(ID).p(1).coord(1) DEV(ID).p(2).coord(1) DEV(ID).p(3).coord(1) DEV(ID).p(4).coord(1) ], ...
                            [DEV(ID).p(1).coord(2) DEV(ID).p(2).coord(2) DEV(ID).p(3).coord(2) DEV(ID).p(4).coord(2) ])

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
        % mat.n_of_modes = n_of_modes;
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

        case 'TAPSQDGFETLONG'
        % n_of_modes = 1;
        mat.type = 'TAPSQDGFETLONG';
        dx = arg_dx; % transport direction
        dy = arg_dy; % confinement direction

        n_of_valleys = 1;  % number of valleys
        deg_factor = 1;    % degenarcy factor of the valley
        
        L_c = 10; % channel length
        L_s = 30; % source length
        L_d = 30; % drain length
        
        L_x = L_c + L_s + L_d;
        
        W_c  = arg_w_ch; % channel width

        SQU  = arg_offsetU; % squeeze upper
        SQL  = arg_offsetL; % squeeze lower

        taper = arg_taper;

        W_ox = arg_w_ox; % oxid width
        
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
        % me_x_ox = 0.041;
        % me_y_ox = 0.041;
        % me_z_ox = 0.041;
        
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
        DEV(1).p(2).coord = [L_s-taper/2 0];
        DEV(1).p(3).coord = [L_s-taper/2 W_c];
        DEV(1).p(4).coord = [0 W_c];
        
        DEV(2).Eg = Eg_c;
        DEV(2).Eps = Eps_c;
        DEV(2).me_x = me_x_ch;
        DEV(2).me_y = me_y_ch;
        DEV(2).me_z = me_z_ch;
        DEV(2).V = Eg_c;
        DEV(2).Nd = 0;
        DEV(2).p(1).coord = [L_s-taper/2 0+SQL*0];
        DEV(2).p(2).coord = [L_s+taper/2 0+SQL*1];
        DEV(2).p(3).coord = [L_s+taper/2 W_c-SQU*1];
        DEV(2).p(4).coord = [L_s-taper/2 W_c-SQU*0];



        DEV(3).Eg = Eg_c;
        DEV(3).Eps = Eps_c;
        DEV(3).me_x = me_x_ch;
        DEV(3).me_y = me_y_ch;
        DEV(3).me_z = me_z_ch;
        DEV(3).V = Eg_c;
        DEV(3).Nd = 0;
        DEV(3).p(1).coord = [L_s+taper/2 0+SQL];
        DEV(3).p(2).coord = [L_s+L_c-taper/2 0+SQL];
        DEV(3).p(3).coord = [L_s+L_c-taper/2 W_c-SQU];
        DEV(3).p(4).coord = [L_s+taper/2 W_c-SQU];

        DEV(4).Eg = Eg_c;
        DEV(4).Eps = Eps_c;
        DEV(4).me_x = me_x_ch;
        DEV(4).me_y = me_y_ch;
        DEV(4).me_z = me_z_ch;
        DEV(4).V = Eg_c;
        DEV(4).Nd = 0;
        DEV(4).p(1).coord = [L_s+L_c-taper/2 0+SQL*1];
        DEV(4).p(2).coord = [L_s+L_c+taper/2 0+SQL*0];
        DEV(4).p(3).coord = [L_s+L_c+taper/2 W_c-SQU*0];
        DEV(4).p(4).coord = [L_s+L_c-taper/2 W_c-SQU*1];
        
        DEV(5).Eg = Eg_c;
        DEV(5).Eps = Eps_c;
        DEV(5).me_x = me_x_ch;
        DEV(5).me_y = me_y_ch;
        DEV(5).me_z = me_z_ch;
        DEV(5).V = Eg_c;
        DEV(5).Nd = N_d;
        DEV(5).p(1).coord = [L_s+L_c+taper/2 0];
        DEV(5).p(2).coord = [L_s+L_c+L_d 0];
        DEV(5).p(3).coord = [L_s+L_c+L_d W_c];
        DEV(5).p(4).coord = [L_s+L_c+taper/2 W_c];
        
        DEV(6).Eg = Eg_ox;
        DEV(6).Eps = Eps_ox;
        DEV(6).me_x = me_x_ox;
        DEV(6).me_y = me_y_ox;
        DEV(6).me_z = me_z_ox;
        DEV(6).V = Eg_c+dEc*(Eg_ox-Eg_c);
        DEV(6).Nd = 0;
        DEV(6).p(1).coord = [0 -W_ox];
        DEV(6).p(2).coord = [L_x -W_ox];
        DEV(6).p(3).coord = [L_x W_c+W_ox];
        DEV(6).p(4).coord = [0 W_c+W_ox];
        
        % DEV(5).Eg = Eg_ox;
        % DEV(5).Eps = Eps_ox;
        % DEV(5).me_x = me_x_ox;
        % DEV(5).me_y = me_y_ox;
        % DEV(5).me_z = me_z_ox;
        % DEV(5).V = Eg_c+dEc*(Eg_ox-Eg_c);
        % DEV(5).Nd = 0;
        % DEV(5).p(1).coord = [0 W_c];
        % DEV(5).p(2).coord = [L_x W_c];
        % DEV(5).p(3).coord = [L_x W_c+W_ox];
        % DEV(5).p(4).coord = [0 W_c+W_c];
        % 
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
        
        DEVo = [6, 2, 4,3,1,5];
        for IX = 1 : mat.Nx
            for IY = 1 : mat.Ny
                for IDx = 1 : length(DEV)
                    ID = DEVo(IDx);
                    % if mat.x(IX) >= DEV(ID).p(1).coord(1) && mat.x(IX) <= DEV(ID).p(2).coord(1) && ...
                    %         mat.y(IY) >=DEV(ID).p(2).coord(2) && mat.y(IY) <= DEV(ID).p(3).coord(2)
                    if inpolygon(mat.x(IX), mat.y(IY),[DEV(ID).p(1).coord(1) DEV(ID).p(2).coord(1) DEV(ID).p(3).coord(1) DEV(ID).p(4).coord(1) ], ...
                            [DEV(ID).p(1).coord(2) DEV(ID).p(2).coord(2) DEV(ID).p(3).coord(2) DEV(ID).p(4).coord(2) ])

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
        % mat.n_of_modes = n_of_modes;
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


        case 'TAPSQDGFET2'
        % n_of_modes = 1;
        mat.type = 'TAPSQDGFET2';
        dx = arg_dx; % transport direction
        dy = arg_dy; % confinement direction

        n_of_valleys = 1;  % number of valleys
        deg_factor = 1;    % degenarcy factor of the valley
        
        L_c = 10; % channel length
        L_s = 20; % source length
        L_d = 20; % drain length
        
        L_x = L_c + L_s + L_d;
        
        W_c  = arg_w_ch; % channel width

        SQU  = arg_offsetU; % squeeze upper
        SQL  = arg_offsetL; % squeeze lower

        taper = arg_taper;

        W_ox = 1; % oxid width
        
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
        DEV(1).p(2).coord = [L_s-taper/2 0];
        DEV(1).p(3).coord = [L_s-taper/2 W_c];
        DEV(1).p(4).coord = [0 W_c];
        
        DEV(2).Eg = Eg_c;
        DEV(2).Eps = Eps_c;
        DEV(2).me_x = me_x_ch;
        DEV(2).me_y = me_y_ch;
        DEV(2).me_z = me_z_ch;
        DEV(2).V = Eg_c;
        DEV(2).Nd = 0;
        DEV(2).p(1).coord = [L_s-taper/2 0+SQL*0];
        DEV(2).p(2).coord = [L_s+taper/2 0+SQL*1];
        DEV(2).p(3).coord = [L_s+taper/2 W_c-SQU*1];
        DEV(2).p(4).coord = [L_s-taper/2 W_c-SQU*0];



        DEV(3).Eg = Eg_c;
        DEV(3).Eps = Eps_c;
        DEV(3).me_x = me_x_ch;
        DEV(3).me_y = me_y_ch;
        DEV(3).me_z = me_z_ch;
        DEV(3).V = Eg_c;
        DEV(3).Nd = 0;
        DEV(3).p(1).coord = [L_s+taper/2 0+SQL];
        DEV(3).p(2).coord = [L_s+L_c-taper/2 0+SQL];
        DEV(3).p(3).coord = [L_s+L_c-taper/2 W_c-SQU];
        DEV(3).p(4).coord = [L_s+taper/2 W_c-SQU];

        DEV(4).Eg = Eg_c;
        DEV(4).Eps = Eps_c;
        DEV(4).me_x = me_x_ch;
        DEV(4).me_y = me_y_ch;
        DEV(4).me_z = me_z_ch;
        DEV(4).V = Eg_c;
        DEV(4).Nd = 0;
        DEV(4).p(1).coord = [L_s+L_c-taper/2 0+SQL*1];
        DEV(4).p(2).coord = [L_s+L_c+taper/2 0+SQL*0];
        DEV(4).p(3).coord = [L_s+L_c+taper/2 W_c-SQU*0];
        DEV(4).p(4).coord = [L_s+L_c-taper/2 W_c-SQU*1];
        
        DEV(5).Eg = Eg_c;
        DEV(5).Eps = Eps_c;
        DEV(5).me_x = me_x_ch;
        DEV(5).me_y = me_y_ch;
        DEV(5).me_z = me_z_ch;
        DEV(5).V = Eg_c;
        DEV(5).Nd = N_d;
        DEV(5).p(1).coord = [L_s+L_c+taper/2 0];
        DEV(5).p(2).coord = [L_s+L_c+L_d 0];
        DEV(5).p(3).coord = [L_s+L_c+L_d W_c];
        DEV(5).p(4).coord = [L_s+L_c+taper/2 W_c];
        
        DEV(6).Eg = Eg_ox;
        DEV(6).Eps = Eps_ox;
        DEV(6).me_x = me_x_ox;
        DEV(6).me_y = me_y_ox;
        DEV(6).me_z = me_z_ox;
        DEV(6).V = Eg_c+dEc*(Eg_ox-Eg_c);
        DEV(6).Nd = 0;
        DEV(6).p(1).coord = [0 -W_ox];
        DEV(6).p(2).coord = [L_x -W_ox];
        DEV(6).p(3).coord = [L_x W_c+W_ox];
        DEV(6).p(4).coord = [0 W_c+W_ox];
        
        % DEV(5).Eg = Eg_ox;
        % DEV(5).Eps = Eps_ox;
        % DEV(5).me_x = me_x_ox;
        % DEV(5).me_y = me_y_ox;
        % DEV(5).me_z = me_z_ox;
        % DEV(5).V = Eg_c+dEc*(Eg_ox-Eg_c);
        % DEV(5).Nd = 0;
        % DEV(5).p(1).coord = [0 W_c];
        % DEV(5).p(2).coord = [L_x W_c];
        % DEV(5).p(3).coord = [L_x W_c+W_ox];
        % DEV(5).p(4).coord = [0 W_c+W_c];
        % 
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
        
        DEVo = [6, 2, 4,3,1,5];
        for IX = 1 : mat.Nx
            for IY = 1 : mat.Ny
                for IDx = 1 : length(DEV)
                    ID = DEVo(IDx);
                    % if mat.x(IX) >= DEV(ID).p(1).coord(1) && mat.x(IX) <= DEV(ID).p(2).coord(1) && ...
                    %         mat.y(IY) >=DEV(ID).p(2).coord(2) && mat.y(IY) <= DEV(ID).p(3).coord(2)
                    if inpolygon(mat.x(IX), mat.y(IY),[DEV(ID).p(1).coord(1) DEV(ID).p(2).coord(1) DEV(ID).p(3).coord(1) DEV(ID).p(4).coord(1) ], ...
                            [DEV(ID).p(1).coord(2) DEV(ID).p(2).coord(2) DEV(ID).p(3).coord(2) DEV(ID).p(4).coord(2) ])

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
        % mat.n_of_modes = n_of_modes;
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


        case 'TAPSQDGFETFORCOMP'
        mat.type = 'TAPSQDGFETFORCOMP';    
        % n_of_modes = 1;

        dx = 0.25; % transport direction
        dy = 0.25; % confinement direction

        n_of_valleys = 1;  % number of valleys
        deg_factor = 1;    % degenarcy factor of the valley
        
        L_c = 10; % channel length
        L_s = 20; % source length
        L_d = 20; % drain length
        
        L_x = L_c + L_s + L_d;
        
        W_c  = 2; % channel width

        SQU  = 0*0.5; % squeeze upper
        SQL  = 0*0.5; % squeeze lower

        taper = 0*0.5;

        W_ox = 1.5; % oxid width
        
        L_y = W_c + W_ox + W_ox;
        
        N_s = 2E+19; % source doping [cm^-3]
        N_d = 2E+19; % drain doping [cm^-3];
        
        Eg_c = 0.74; % bandgap of InGaAs
        Eps_c = 13.9;
        
        me_x_ch = 0.041;
        me_y_ch = 0.041;
        me_z_ch = 0.041;
        
        Eg_ox = 1e6;
        Eps_ox = 3.9;
        
        me_x_ox = 0.041;
        me_y_ox = 0.041;
        me_z_ox = 0.041;
        
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
        DEV(1).p(2).coord = [L_s-taper/2 0];
        DEV(1).p(3).coord = [L_s-taper/2 W_c];
        DEV(1).p(4).coord = [0 W_c];
        
        DEV(2).Eg = Eg_c;
        DEV(2).Eps = Eps_c;
        DEV(2).me_x = me_x_ch;
        DEV(2).me_y = me_y_ch;
        DEV(2).me_z = me_z_ch;
        DEV(2).V = Eg_c;
        DEV(2).Nd = 0;
        DEV(2).p(1).coord = [L_s-taper/2 0+SQL*0];
        DEV(2).p(2).coord = [L_s+taper/2 0+SQL*1];
        DEV(2).p(3).coord = [L_s+taper/2 W_c-SQU*1];
        DEV(2).p(4).coord = [L_s-taper/2 W_c-SQU*0];



        DEV(3).Eg = Eg_c;
        DEV(3).Eps = Eps_c;
        DEV(3).me_x = me_x_ch;
        DEV(3).me_y = me_y_ch;
        DEV(3).me_z = me_z_ch;
        DEV(3).V = Eg_c;
        DEV(3).Nd = 0;
        DEV(3).p(1).coord = [L_s+taper/2 0+SQL];
        DEV(3).p(2).coord = [L_s+L_c-taper/2 0+SQL];
        DEV(3).p(3).coord = [L_s+L_c-taper/2 W_c-SQU];
        DEV(3).p(4).coord = [L_s+taper/2 W_c-SQU];

        DEV(4).Eg = Eg_c;
        DEV(4).Eps = Eps_c;
        DEV(4).me_x = me_x_ch;
        DEV(4).me_y = me_y_ch;
        DEV(4).me_z = me_z_ch;
        DEV(4).V = Eg_c;
        DEV(4).Nd = 0;
        DEV(4).p(1).coord = [L_s+L_c-taper/2 0+SQL*1];
        DEV(4).p(2).coord = [L_s+L_c+taper/2 0+SQL*0];
        DEV(4).p(3).coord = [L_s+L_c+taper/2 W_c-SQU*0];
        DEV(4).p(4).coord = [L_s+L_c-taper/2 W_c-SQU*1];
        
        DEV(5).Eg = Eg_c;
        DEV(5).Eps = Eps_c;
        DEV(5).me_x = me_x_ch;
        DEV(5).me_y = me_y_ch;
        DEV(5).me_z = me_z_ch;
        DEV(5).V = Eg_c;
        DEV(5).Nd = N_d;
        DEV(5).p(1).coord = [L_s+L_c+taper/2 0];
        DEV(5).p(2).coord = [L_s+L_c+L_d 0];
        DEV(5).p(3).coord = [L_s+L_c+L_d W_c];
        DEV(5).p(4).coord = [L_s+L_c+taper/2 W_c];
        
        DEV(6).Eg = Eg_ox;
        DEV(6).Eps = Eps_ox;
        DEV(6).me_x = me_x_ox;
        DEV(6).me_y = me_y_ox;
        DEV(6).me_z = me_z_ox;
        DEV(6).V = Eg_c+dEc*(Eg_ox-Eg_c);
        DEV(6).Nd = 0;
        DEV(6).p(1).coord = [0 -W_ox];
        DEV(6).p(2).coord = [L_x -W_ox];
        DEV(6).p(3).coord = [L_x W_c+W_ox];
        DEV(6).p(4).coord = [0 W_c+W_ox];
        
        % DEV(5).Eg = Eg_ox;
        % DEV(5).Eps = Eps_ox;
        % DEV(5).me_x = me_x_ox;
        % DEV(5).me_y = me_y_ox;
        % DEV(5).me_z = me_z_ox;
        % DEV(5).V = Eg_c+dEc*(Eg_ox-Eg_c);
        % DEV(5).Nd = 0;
        % DEV(5).p(1).coord = [0 W_c];
        % DEV(5).p(2).coord = [L_x W_c];
        % DEV(5).p(3).coord = [L_x W_c+W_ox];
        % DEV(5).p(4).coord = [0 W_c+W_c];
        % 
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
        
        DEVo = [6, 2, 4,3,1,5];
        for IX = 1 : mat.Nx
            for IY = 1 : mat.Ny
                for IDx = 1 : length(DEV)
                    ID = DEVo(IDx);
                    % if mat.x(IX) >= DEV(ID).p(1).coord(1) && mat.x(IX) <= DEV(ID).p(2).coord(1) && ...
                    %         mat.y(IY) >=DEV(ID).p(2).coord(2) && mat.y(IY) <= DEV(ID).p(3).coord(2)
                    if inpolygon(mat.x(IX), mat.y(IY),[DEV(ID).p(1).coord(1) DEV(ID).p(2).coord(1) DEV(ID).p(3).coord(1) DEV(ID).p(4).coord(1) ], ...
                            [DEV(ID).p(1).coord(2) DEV(ID).p(2).coord(2) DEV(ID).p(3).coord(2) DEV(ID).p(4).coord(2) ])

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
        % mat.n_of_modes = n_of_modes;
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

        case 'TAPSQDGFETFORCOMPMV'
            mat.type = 'TAPSQDGFETFORCOMPMV';    
        % n_of_modes = 1;

        dx = arg_dx; % transport direction
        dy = arg_dy; % confinement direction

        n_of_valleys = 1;  % number of valleys
        deg_factor = 1;    % degenarcy factor of the valley
        
        L_c = 10; % channel length
        L_s = 20; % source length
        L_d = 20; % drain length
        
        L_x = L_c + L_s + L_d;
        
        W_c  = arg_w_ch; % channel width

        SQU  = 0*0.5; % squeeze upper
        SQL  = 0*0.5; % squeeze lower

        taper = 0*0.5;

        W_ox = arg_w_ox; % oxid width
        
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
        DEV(1).p(2).coord = [L_s-taper/2 0];
        DEV(1).p(3).coord = [L_s-taper/2 W_c];
        DEV(1).p(4).coord = [0 W_c];
        
        DEV(2).Eg = Eg_c;
        DEV(2).Eps = Eps_c;
        DEV(2).me_x = me_x_ch;
        DEV(2).me_y = me_y_ch;
        DEV(2).me_z = me_z_ch;
        DEV(2).V = Eg_c;
        DEV(2).Nd = 0;
        DEV(2).p(1).coord = [L_s-taper/2 0+SQL*0];
        DEV(2).p(2).coord = [L_s+taper/2 0+SQL*1];
        DEV(2).p(3).coord = [L_s+taper/2 W_c-SQU*1];
        DEV(2).p(4).coord = [L_s-taper/2 W_c-SQU*0];



        DEV(3).Eg = Eg_c;
        DEV(3).Eps = Eps_c;
        DEV(3).me_x = me_x_ch;
        DEV(3).me_y = me_y_ch;
        DEV(3).me_z = me_z_ch;
        DEV(3).V = Eg_c;
        DEV(3).Nd = 0;
        DEV(3).p(1).coord = [L_s+taper/2 0+SQL];
        DEV(3).p(2).coord = [L_s+L_c-taper/2 0+SQL];
        DEV(3).p(3).coord = [L_s+L_c-taper/2 W_c-SQU];
        DEV(3).p(4).coord = [L_s+taper/2 W_c-SQU];

        DEV(4).Eg = Eg_c;
        DEV(4).Eps = Eps_c;
        DEV(4).me_x = me_x_ch;
        DEV(4).me_y = me_y_ch;
        DEV(4).me_z = me_z_ch;
        DEV(4).V = Eg_c;
        DEV(4).Nd = 0;
        DEV(4).p(1).coord = [L_s+L_c-taper/2 0+SQL*1];
        DEV(4).p(2).coord = [L_s+L_c+taper/2 0+SQL*0];
        DEV(4).p(3).coord = [L_s+L_c+taper/2 W_c-SQU*0];
        DEV(4).p(4).coord = [L_s+L_c-taper/2 W_c-SQU*1];
        
        DEV(5).Eg = Eg_c;
        DEV(5).Eps = Eps_c;
        DEV(5).me_x = me_x_ch;
        DEV(5).me_y = me_y_ch;
        DEV(5).me_z = me_z_ch;
        DEV(5).V = Eg_c;
        DEV(5).Nd = N_d;
        DEV(5).p(1).coord = [L_s+L_c+taper/2 0];
        DEV(5).p(2).coord = [L_s+L_c+L_d 0];
        DEV(5).p(3).coord = [L_s+L_c+L_d W_c];
        DEV(5).p(4).coord = [L_s+L_c+taper/2 W_c];
        
        DEV(6).Eg = Eg_ox;
        DEV(6).Eps = Eps_ox;
        DEV(6).me_x = me_x_ox;
        DEV(6).me_y = me_y_ox;
        DEV(6).me_z = me_z_ox;
        DEV(6).V = Eg_c+dEc*(Eg_ox-Eg_c);
        DEV(6).Nd = 0;
        DEV(6).p(1).coord = [0 -W_ox];
        DEV(6).p(2).coord = [L_x -W_ox];
        DEV(6).p(3).coord = [L_x W_c+W_ox];
        DEV(6).p(4).coord = [0 W_c+W_ox];
        
        % DEV(5).Eg = Eg_ox;
        % DEV(5).Eps = Eps_ox;
        % DEV(5).me_x = me_x_ox;
        % DEV(5).me_y = me_y_ox;
        % DEV(5).me_z = me_z_ox;
        % DEV(5).V = Eg_c+dEc*(Eg_ox-Eg_c);
        % DEV(5).Nd = 0;
        % DEV(5).p(1).coord = [0 W_c];
        % DEV(5).p(2).coord = [L_x W_c];
        % DEV(5).p(3).coord = [L_x W_c+W_ox];
        % DEV(5).p(4).coord = [0 W_c+W_c];
        % 
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
        
        DEVo = [6, 2, 4,3,1,5];
        for IX = 1 : mat.Nx
            for IY = 1 : mat.Ny
                for IDx = 1 : length(DEV)
                    ID = DEVo(IDx);
                    % if mat.x(IX) >= DEV(ID).p(1).coord(1) && mat.x(IX) <= DEV(ID).p(2).coord(1) && ...
                    %         mat.y(IY) >=DEV(ID).p(2).coord(2) && mat.y(IY) <= DEV(ID).p(3).coord(2)
                    if inpolygon(mat.x(IX), mat.y(IY),[DEV(ID).p(1).coord(1) DEV(ID).p(2).coord(1) DEV(ID).p(3).coord(1) DEV(ID).p(4).coord(1) ], ...
                            [DEV(ID).p(1).coord(2) DEV(ID).p(2).coord(2) DEV(ID).p(3).coord(2) DEV(ID).p(4).coord(2) ])

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
        % mat.n_of_modes = n_of_modes;
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
        
        case 'TAPSQDGFETFORCOMPMV2'
            mat.type = 'TAPSQDGFETFORCOMPMV2';    
        % n_of_modes = 1;

        dx = arg_dx; % transport direction
        dy = arg_dy; % confinement direction

        n_of_valleys = 1;  % number of valleys
        deg_factor = 1;    % degenarcy factor of the valley
        
        L_c = 10; % channel length
        L_s = 20; % source length
        L_d = 20; % drain length
        
        L_x = L_c + L_s + L_d;
        
        W_c  = arg_w_ch; % channel width

        SQU  = 0.75; % squeeze upper
        SQL  = 0.75; % squeeze lower

        taper = 0.5;

        W_ox = arg_w_ox; % oxid width
        
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
        DEV(1).p(2).coord = [L_s-taper/2 0];
        DEV(1).p(3).coord = [L_s-taper/2 W_c];
        DEV(1).p(4).coord = [0 W_c];
        
        DEV(2).Eg = Eg_c;
        DEV(2).Eps = Eps_c;
        DEV(2).me_x = me_x_ch;
        DEV(2).me_y = me_y_ch;
        DEV(2).me_z = me_z_ch;
        DEV(2).V = Eg_c;
        DEV(2).Nd = 0;
        DEV(2).p(1).coord = [L_s-taper/2 0+SQL*0];
        DEV(2).p(2).coord = [L_s+taper/2 0+SQL*1];
        DEV(2).p(3).coord = [L_s+taper/2 W_c-SQU*1];
        DEV(2).p(4).coord = [L_s-taper/2 W_c-SQU*0];



        DEV(3).Eg = Eg_c;
        DEV(3).Eps = Eps_c;
        DEV(3).me_x = me_x_ch;
        DEV(3).me_y = me_y_ch;
        DEV(3).me_z = me_z_ch;
        DEV(3).V = Eg_c;
        DEV(3).Nd = 0;
        DEV(3).p(1).coord = [L_s+taper/2 0+SQL];
        DEV(3).p(2).coord = [L_s+L_c-taper/2 0+SQL];
        DEV(3).p(3).coord = [L_s+L_c-taper/2 W_c-SQU];
        DEV(3).p(4).coord = [L_s+taper/2 W_c-SQU];

        DEV(4).Eg = Eg_c;
        DEV(4).Eps = Eps_c;
        DEV(4).me_x = me_x_ch;
        DEV(4).me_y = me_y_ch;
        DEV(4).me_z = me_z_ch;
        DEV(4).V = Eg_c;
        DEV(4).Nd = 0;
        DEV(4).p(1).coord = [L_s+L_c-taper/2 0+SQL*1];
        DEV(4).p(2).coord = [L_s+L_c+taper/2 0+SQL*0];
        DEV(4).p(3).coord = [L_s+L_c+taper/2 W_c-SQU*0];
        DEV(4).p(4).coord = [L_s+L_c-taper/2 W_c-SQU*1];
        
        DEV(5).Eg = Eg_c;
        DEV(5).Eps = Eps_c;
        DEV(5).me_x = me_x_ch;
        DEV(5).me_y = me_y_ch;
        DEV(5).me_z = me_z_ch;
        DEV(5).V = Eg_c;
        DEV(5).Nd = N_d;
        DEV(5).p(1).coord = [L_s+L_c+taper/2 0];
        DEV(5).p(2).coord = [L_s+L_c+L_d 0];
        DEV(5).p(3).coord = [L_s+L_c+L_d W_c];
        DEV(5).p(4).coord = [L_s+L_c+taper/2 W_c];
        
        DEV(6).Eg = Eg_ox;
        DEV(6).Eps = Eps_ox;
        DEV(6).me_x = me_x_ox;
        DEV(6).me_y = me_y_ox;
        DEV(6).me_z = me_z_ox;
        DEV(6).V = Eg_c+dEc*(Eg_ox-Eg_c);
        DEV(6).Nd = 0;
        DEV(6).p(1).coord = [0 -W_ox];
        DEV(6).p(2).coord = [L_x -W_ox];
        DEV(6).p(3).coord = [L_x W_c+W_ox];
        DEV(6).p(4).coord = [0 W_c+W_ox];
        
        % DEV(5).Eg = Eg_ox;
        % DEV(5).Eps = Eps_ox;
        % DEV(5).me_x = me_x_ox;
        % DEV(5).me_y = me_y_ox;
        % DEV(5).me_z = me_z_ox;
        % DEV(5).V = Eg_c+dEc*(Eg_ox-Eg_c);
        % DEV(5).Nd = 0;
        % DEV(5).p(1).coord = [0 W_c];
        % DEV(5).p(2).coord = [L_x W_c];
        % DEV(5).p(3).coord = [L_x W_c+W_ox];
        % DEV(5).p(4).coord = [0 W_c+W_c];
        % 
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
        
        DEVo = [6, 2, 4,3,1,5];
        for IX = 1 : mat.Nx
            for IY = 1 : mat.Ny
                for IDx = 1 : length(DEV)
                    ID = DEVo(IDx);
                    % if mat.x(IX) >= DEV(ID).p(1).coord(1) && mat.x(IX) <= DEV(ID).p(2).coord(1) && ...
                    %         mat.y(IY) >=DEV(ID).p(2).coord(2) && mat.y(IY) <= DEV(ID).p(3).coord(2)
                    if inpolygon(mat.x(IX), mat.y(IY),[DEV(ID).p(1).coord(1) DEV(ID).p(2).coord(1) DEV(ID).p(3).coord(1) DEV(ID).p(4).coord(1) ], ...
                            [DEV(ID).p(1).coord(2) DEV(ID).p(2).coord(2) DEV(ID).p(3).coord(2) DEV(ID).p(4).coord(2) ])

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
        % mat.n_of_modes = n_of_modes;
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
end


