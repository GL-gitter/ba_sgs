function dc_dEf = eval_slice_density_derivative(Esub, Ef,mat)

h = 6.626E-34;
hb = h/2/pi;
m0 = 9.109E-31;
kB = 1.38E-23;
q = 1.602E-19;

Temp = mat.Temp;

me_x = mat.me_x_ch;
me_z = mat.me_z_ch;

n_of_valleys = mat.n_of_valleys;
n_of_modes = mat.n_of_modes;
deg_factor = mat.deg_factor;

% effective density of states mass

m_dos = (me_x.*me_z)^(1/2);

NC = m0*m_dos*kB*Temp/(pi*hb.^2);

dc_dEf = 0;

for IV = 1 : n_of_valleys
    for IM = 1 : n_of_modes
        
        dc_dEf_temp = deg_factor(IV).*NC(IV).*(q/(kB*Temp))*...
            (1./(1+exp(q*(Esub(IV,IM)-Ef)/(kB*Temp))));
        
        dc_dEf = dc_dEf + dc_dEf_temp;
        
    end
end

