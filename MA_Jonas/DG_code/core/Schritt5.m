% Erst "Schritt3.m" und
% "Schritt4.m" für die nötigen Parameter in dem Workspace ausführen
%% Definitionen
sigma_x = [0 1; 1 0]; % Pauli-Matrix
sigma_y = [0 -1i; 1i 0];
sigma_z = [1 0; 0 -1];
mu = 0.0; % eV, =0 für Dirac-Punkt
Delta0 = 0.002; % eV
L = 10e-9; % 10 nm für L<<xi_0
U0 = 0; % eV
k_y_werte = linspace(-3/L,3/L,31); % 1/L
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
        if abs(sin(theta_e(j))) < bed_winkel(E)
            theta_e_gueltig(i,j) = theta_e(j);
        end
    end
end

% ABS berechnen, numerisch
% M_L = @(E,phi) exp(1i*phi/2)*exp(1i*beta*sigma_x);
% M_R = @(E,phi) exp(-1i*phi/2)*exp(-1i*beta*sigma_x);
det_ABS = @(E,phi,k_y) det(T_h(E,k_y)*M_L(E,phi) - M_R(E,phi)*T_e(E,k_y));
abs_det_ABS = @(E,phi,k_y) abs(det_ABS(E,phi,k_y));
% Determinante
% det_werte = zeros(length(E_werte), length(phi_werte), length(k_y_werte));
% for m = 1:length(k_y_werte)
%     k_y = k_y_werte(m);
%     for j = 1:length(phi_werte)
%         phi = phi_werte(j);
%         for i = 1:length(E_werte)
%             E = E_werte(i);
%             det_werte(i,j,m) = abs_det_ABS(E,phi,k_y);
%         end
%     end
% end

E_ABS_plus_num = zeros(length(phi_werte),length(k_y_werte));
E_ABS_minus_num = zeros(length(phi_werte),length(k_y_werte));
options = optimset('TolX', 1e-12);
for m = 1:length(k_y_werte)
    k_y = k_y_werte(m);
    for j = 1:length(phi_werte)
        phi = phi_werte(j);
        % Lösen von det(...)=0
        E_ABS_plus_num(j,m) = fminbnd(@(E) abs_det_ABS(E,phi,k_y), 0, Delta0, options);
        E_ABS_minus_num(j,m) = fminbnd(@(E) abs_det_ABS(E,phi,k_y), -Delta0, 0, options);
    end
end
% ABS berechnen, analytisch
E_ABS_plus_ana = zeros(length(phi_werte),length(k_y_werte));
E_ABS_minus_ana = zeros(length(phi_werte),length(k_y_werte));
DP_ana = @(phi, k_y) Delta0*sqrt(1-sech(k_y*L)^2*sin(phi/2)^2);
for m = 1:length(k_y_werte)
    k_y = k_y_werte(m);
    for j = 1:length(phi_werte)
        phi = phi_werte(j);
        E_ABS_plus_ana(j,m) = DP_ana(phi,k_y);
        E_ABS_minus_ana(j,m) = -DP_ana(phi,k_y);
    end
end


% Darstellung der Energien
% Analytische Lösungen
figure
mesh(k_y_werte, phi_werte, E_ABS_plus_ana)
hold on
mesh(k_y_werte, phi_werte, E_ABS_minus_ana)
xlabel('k_y (m^{-1})')
ylabel('\phi (rad)')
zlabel('E_{ABS} (eV)')
% legend('Numerische Lösung E_+', 'Numerische Lösung E_-', ...
%     'Analytische Lösung E_+', 'Analytische Lösung E_-')
grid on
title('ABS-Analytische Lösungen')
% Numerische Lösungen
figure
mesh(k_y_werte, phi_werte, E_ABS_plus_num)
hold on
mesh(k_y_werte, phi_werte, E_ABS_minus_num)
xlabel('k_y (m^{-1})')
ylabel('\phi (rad)')
zlabel('E_{ABS} (eV)')
% legend('Numerische Lösung E_+', 'Numerische Lösung E_-')
grid on
title('ABS-Numerische Lösungen')
% Fehler von numerischer zu analytischer Rechnung (Für höhere Präzision 
% kleinere L wählen)
Fehler_ABS_plus = abs(E_ABS_plus_num - E_ABS_plus_ana);
Fehler_ABS_minus = abs(E_ABS_minus_num - E_ABS_minus_ana);

figure
mesh(k_y_werte, phi_werte, Fehler_ABS_plus)
xlabel('k_y (m^{-1})')
ylabel('\phi (rad)')
zlabel('Absoluter Fehler (eV)')
title('Abweichung von der analytischen ABS-Lösung')
grid on

% SAR und RAR verifizieren
mu_werte = [0.1, 0]; % RAR und SAR
v_y_h_num = NaN(length(E_werte),length(theta_e),length(mu_werte));
v_y_h_ana = NaN(length(E_werte),length(theta_e),length(mu_werte));
for j = 1:length(E_werte)
    E_verif = E_werte(j);
    if E_verif <= 0
        continue
    end
    for l=1:length(theta_e)
        theta_verif = theta_e(l);
        for i = 1:length(mu_werte)
            mu_verif = mu_werte(i);
            % Winkelbedingung überprüfen
            if abs(sin(theta_verif)) >= ...
                    abs(mu_verif-U0-E_verif)/abs(mu_verif-U0+E_verif)
                continue
            end
            % Transversale Wellenzahl
            k_y_verif = abs(mu_verif-U0+E_verif)/(hbar*vF)*sin(theta_verif);
            % Longitudinale Wellenzahl des Lochs
            k_x_verif = sqrt(((mu_verif-U0-E_verif)/(hbar*vF))^2-k_y_verif^2);
            % Lochpropagation
            G_h = k_y_verif*sigma_z + 1i*(mu_verif-U0-E_verif)/(hbar*vF)*sigma_x;
            % Eigenvektoren der Lochmoden
            [EV,EW] = eig(G_h);
            for n = 1:2
                psi_h = EV(:,n);
                % Geschwindigkeiten numerisch bestimmen
                v_x = real(-vF*(psi_h'*sigma_x*psi_h)/(psi_h'*psi_h));
                v_y = real(-vF*(psi_h'*sigma_y*psi_h)/(psi_h'*psi_h));
                if v_x < 0
                    v_y_h_num(j,l,i) = v_y;
                end
            end
            % Analytische Vergleichslösung
            v_y_h_ana(j,l,i) = -sign(mu_verif-U0-E_verif) ...
                *vF*k_y_verif/sqrt(k_x_verif^2+k_y_verif^2);
        end
    end
end
figure
% RAR
subplot(1,2,1)
mesh(theta_e,E_werte,v_y_h_num(:,:,1))
xlabel('\theta_e (rad)')
ylabel('E (eV)')
zlabel('v_{y,h} (m/s)')
title('RAR')
% SAR
subplot(1,2,2)
mesh(theta_e,E_werte,v_y_h_num(:,:,2))
xlabel('\theta_e (rad)')
ylabel('E (eV)')
zlabel('v_{y,h} (m/s)')
title('SAR')