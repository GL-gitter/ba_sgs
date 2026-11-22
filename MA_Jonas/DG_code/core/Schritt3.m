%% Definitionen
sigma_x = [0 1; 1 0]; % Pauli-Matrix
hbar = 6.582119569e-16; % eV*s
mu = 0.1; % eV
vF = 1e6; % m/s
Delta0 = 0.002; % eV
L = 10e-9; % 10 nm für L<<xi_0  
U0 = 0; % eV
U = @(x) U0 + 0.*x; % zunächst konstantes Potential

k_e = @(E) (mu - U0 + E)/(hbar*vF);
k_h = @(E) (mu - U0 - E)/(hbar*vF);
T_e = @(E) (cos(k_e(E)*L)*eye(2)+1i*sin(k_e(E)*L)*sigma_x);
T_h = @(E) (cos(k_h(E)*L)*eye(2)+1i*sin(k_h(E)*L)*sigma_x);
beta = @(E) acos(E/Delta0);
M_L = @(E,phi) exp(1i*phi/2)* (cos(beta(E))*eye(2)+1i*sin(beta(E))*sigma_x);
M_R = @(E,phi) exp(-1i*phi/2)* (cos(-beta(E))*eye(2)+1i*sin(-beta(E))*sigma_x);

%% Berechnen der ABS
det_ABS = @(E,phi) det(T_h(E)*M_L(E,phi) - M_R(E,phi)*T_e(E));
abs_det_ABS = @(E,phi) abs(det_ABS(E,phi));

E_werte = linspace(-Delta0, Delta0, 1000);
det_werte = zeros(size(E_werte));
phi_werte = linspace(0, 2*pi, 701);

for j = 1:length(phi_werte)
    phi = phi_werte(j);
    for i = 1:length(E_werte)
        E = E_werte(i);
        det_werte(i,j) = abs_det_ABS(E,phi);
    end
end
% Darstellung der Determinante
% figure
% plot(E_werte, det_werte(:,5))
% xlabel('E (eV)')
% ylabel('|det|')
% grid on
% Berechnung der phasenabhängigen Energien
E_ABS_plus_num = zeros(size(phi_werte));
E_ABS_minus_num = zeros(size(phi_werte));
options = optimset('TolX', 1e-12);
for j = 1:length(phi_werte)
    phi = phi_werte(j);
    % Lösen von det(...)=0
    E_ABS_plus_num(j) = fminbnd(@(E) abs_det_ABS(E,phi), 0, Delta0, options);
    E_ABS_minus_num(j) = fminbnd(@(E) abs_det_ABS(E,phi), -Delta0, 0, options);
end
E_ABS_plus_ana = zeros(size(phi_werte));
E_ABS_minus_ana = zeros(size(phi_werte));
for j = 1:length(phi_werte)
    phi = phi_werte(j);
    E_ABS_plus_ana(j) = Delta0*sqrt(1-sin(phi/2)^2);
    E_ABS_minus_ana(j) = -Delta0*sqrt(1-sin(phi/2)^2);
end
% Darstellung der Energien
% Vergleich
figure
plot(phi_werte, E_ABS_plus_num, 'Linewidth', 1.5); hold on;
plot(phi_werte, E_ABS_minus_num, 'Linewidth', 1.5); hold on;
plot(phi_werte, E_ABS_plus_ana, '--', 'Linewidth', 1.5); hold on;
plot(phi_werte, E_ABS_minus_ana, '--', 'Linewidth', 1.5);
xlabel('\phi (rad)')
ylabel('E_{ABS} (eV)')
legend('Numerische Lösung E_+', 'Numerische Lösung E_-', ...
'Analytische Lösung E_+', 'Analytische Lösung E_-')
grid on
title('ABS-Vergleich')
% Nur numerische Ergebnisse
figure
plot(phi_werte, E_ABS_plus_num, 'Linewidth', 1.5); hold on;
plot(phi_werte, E_ABS_minus_num, 'Linewidth', 1.5); hold on;
xlabel('\phi (rad)')
ylabel('E_{ABS} (eV)')
legend('Numerische Lösung E_+', 'Numerische Lösung E_-')
grid on
title('ABS-Numerische Lösungen')