% Erst "Schritt3.m" und
% "Schritt4.m" für die nötigen Parameter in dem Workspace ausführen
%% Definitionen
sigma_x = [0 1; 1 0]; % Pauli-Matrix
sigma_z = [1 0; 0 -1];
hbar = 6.582119569e-16; % eV*s
mu = 0.0; % eV, =0 für Dirac-Punkt
vF = 1e6; % m/s
Delta0 = 0.002; % eV
L = 10e-9; % 10 nm für L<<xi_0
U0 = 0; % eV
% Gültiges k_y berechnen ???
k_y = @(E) abs((mu-U0+E)/(hbar*vF))*sin(theta_e);
% Wellenvektoren in x-Richtung
k_x_e = @(E,k_y) sqrt(((mu-U0+E)/(hbar*vF))^2 - k_y^2);
k_x_h = @(E,k_y) sqrt(((mu-U0-E)/(hbar*vF))^2 - k_y^2);

% Propagationsmatrizen und Winkelabhängigkeiten
T_e = @(E,k_y) expm(L*(k_y*sigma_z+1i*(E+mu-U0)/(hbar*vF) * sigma_x));
T_h = @(E,k_y) expm(L*(k_y*sigma_z+1i*(mu-U0-E)/(hbar*vF) * sigma_x));
bed_winkel = @(E) abs(mu-U0-E)/abs(mu-U0+E);
theta_e = linspace(-pi/2, pi/2, 1000);
theta_e_gueltig = NaN(length(E_werte),length(theta_e)); % NaN, da 0 auch ein 
                                                        % gültiger Winkel wäre
for i=1:length(E_werte)
    E = E_werte(i);
    for j=1:length(theta_e)
        if abs(sin(theta_e(i))) < bed_winkel(E)
            theta_e_gueltig(i,j) = theta_e(j);
        end
    end
end

% ABS berechnen, numerisch
% M_L = @(E,phi) exp(1i*phi/2)*exp(1i*beta*sigma_x);
% M_R = @(E,phi) exp(-1i*phi/2)*exp(-1i*beta*sigma_x);
det_ABS = @(E,phi,k_y) det(T_h(E,k_y)*M_L(E,phi) - M_R(E,phi)*T_e(E,k_y));
abs_det_ABS = @(E,phi,k_y) abs(det_ABS(E,phi,k_y));
for j = 1:length(phi_werte)
    phi = phi_werte(j);
    for i = 1:length(E_werte)
        E = E_werte(i);
        det_werte(i,j) = abs_det_ABS(E,phi,k_y);
    end
end
E_ABS_plus_num = zeros(size(phi_werte));
E_ABS_minus_num = zeros(size(phi_werte));
options = optimset('TolX', 1e-12);
for j = 1:length(phi_werte)
    phi = phi_werte(j);
    % Lösen von det(...)=0
    E_ABS_plus_num(j) = fminbnd(@(E) abs_det_ABS(E,phi,k_y), 0, Delta0, options);
    E_ABS_minus_num(j) = fminbnd(@(E) abs_det_ABS(E,phi,k_y), -Delta0, 0, options);
end
% ABS berechnen, analytisch
E_ABS_plus_ana = zeros(size(phi_werte));
E_ABS_minus_ana = zeros(size(phi_werte));
DP_ana = @(phi, k_y) Delta0*sqrt(1-sech(k_y*L)^2*sin(phi/2)^2);
for j = 1:length(phi_werte)
    phi = phi_werte(j);
    E_ABS_plus_ana(j) = DP_ana(phi,k_y);
    E_ABS_minus_ana(j) = -DP_ana(phi,k_y);
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

% SAR und RAR verifizieren