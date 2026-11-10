function H_D = eval_Dirac_Hamiltonian(k_x, v_F)
% Nur Betrachtung des K-Valleys (k_y = 0)
% H_D = hb v_F [0, 1; 1, 0] und in Impulsdarstellung mit Pauli-Matrix sigma_x
% H_D = -i hb v_F sigma_x del_x

h = 6.626E-34;
hb = h/2/pi;

sigma_x = [0, 1; 1, 0];

% Mit ebener Welle Psi=ue^(i k_x x):
H_D = hb * v_F * k_x * sigma_x;

end
