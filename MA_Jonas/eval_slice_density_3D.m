function c0 = eval_slice_density_3D(Esub, Ef, mat)

%% Ändern

h = 6.626E-34;
hb = h/2/pi;
q = 1.602E-19;
m0 = 9.109E-31;
kB = 1.38E-23;

Temp = mat.Temp;
me_x = mat.me_x_ch;
me_y = mat.me_y_ch;
me_z = mat.me_z_ch;

n_of_valleys = mat.n_of_valleys;
n_of_modes = mat.n_of_modes;
deg_factor = mat.deg_factor;

% effective density of states

m_dos = (me_x.*me_z)^(1/2);
% NC    = m0*m_dos*kB*Temp/(pi*hb.^2);

c0_mode_valley = zeros(n_of_valleys, n_of_modes);
c0 = 0;

for IV = 1 : n_of_valleys
    for IM = 1 : n_of_modes
        Esubs = squeeze(Esub(IV,IM));
        
        fun = @(x) 2*deg_factor/pi*1./(1+exp(q*((hb*hb*x.*x/2/m0/me_z/q)+Esubs-Ef)/(kB*Temp)));
        % c0_mode_valley(IV, IM) = deg_factor(IV).*NC(IV).*log(1+exp(q*(Ef-Esub(IV,IM))/(kB*Temp)));
        % c0 = c0 + c0_mode_valley(IV, IM);
        c0 = c0 + integral(fun,0,1E10,"RelTol",1e-12);
    end
end
        

