clear; clc;

h = 6.626E-34;
hb = h/2/pi;
v_F = 1e6; % m/s
k_x = linspace(-1e9, 1e9, 201);
eV = 1.602176634e-19;

% Eigenwerte deklarieren
E = zeros(2, length(k_x));
% EW berechnen
for i = 1:length(k_x)
    H_D = eval_Dirac_Hamiltonian(k_x(i),v_F);
    E(:, i) = sort(eig(H_D))/eV;
end

% Analytische Eigenwerte
E_anal_pos = hb * v_F * abs(k_x)/eV; 
E_anal_neg = -hb * v_F * abs(k_x)/eV;

% Plot der beiden Energiebänder
figure;
plot(k_x, E, 'LineWidth', 1.5);
hold on;
plot(k_x, E_anal_pos, 'k--', 'LineWidth', 1.5);
plot(k_x, E_anal_neg, 'k--', 'LineWidth', 1.5);
xlabel('k_x (m^{-1})');
ylabel('Energy (eV)');
grid on;
legend('E_1', 'E_2', 'E_{analytisch,+}', 'E_{analytisch,-}', 'Location', 'best');