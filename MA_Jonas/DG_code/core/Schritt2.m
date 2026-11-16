%% Definitionen
sigma_x = [0 1; 1 0]; % Pauli-Matrix
chi_plus  = [1; 1]/sqrt(2); % Eigenvektoren
chi_minus = [1; -1]/sqrt(2);
% Parameter
hbar = 6.582119569e-16; % eV*s
vF   = 1e6; % m/s
mu = 0.1; % eV
L = 100e-9; % nm
x = linspace(0,L,500);
Delta0 = 0.002; % eV
Phi_S = 0; % Phase
U0 = 0; % eV
U = @(x) U0 + 0.*x; % zunächst konstantes Potential
E_werte  = linspace(-0.99*Delta0, 0.99*Delta0, 200); % eV
r_A_werte = zeros(size(E_werte));
r_N_werte = zeros(size(E_werte));
r_A_ana_werte = zeros(size(E_werte));
A_BdG = vF * blkdiag(sigma_x, -sigma_x);
j_qp = zeros(length(E_werte),length(x));
rho_E_werte = complex(zeros(4,4,length(E_werte)));

% Schleife über die Energien
for i = 1:length(E_werte)
    E_i = E_werte(i);
    beta_i = acos(E_i/Delta0);
    % Wellenzahlen (1/m)
    ke_plus = (mu - U0 + E_i)/(hbar*vF);
    ke_minus = -(mu - U0+ E_i)/(hbar*vF);
    kh_plus = (mu - U0- E_i)/(hbar*vF);
    kh_minus = -(mu - U0- E_i)/(hbar*vF);
    % Moden
    psi_e_plus = chi_plus * exp(1i*ke_plus.*x);
    psi_e_minus = chi_minus * exp(1i*ke_minus.*x);
    psi_h_plus = chi_plus * exp(1i*kh_plus.*x);
    psi_h_minus = chi_minus * exp(1i*kh_minus.*x);

    % Gleichungssystem für die Koeffizienten
    M_S_i = exp(-1i*Phi_S) * (cos(beta_i)*eye(2) - 1i*sin(beta_i)*sigma_x);
    
    A_i = [chi_plus, -M_S_i*chi_minus];
    B_i = M_S_i*chi_plus;
    r_i = A_i\B_i;

    r_A_werte(i) = r_i(1);
    r_N_werte(i) = r_i(2);
    r_A_ana_werte(i) = exp(-1i*(Phi_S + beta_i)); % Bei der Annahme von r_N=0

    % Quasiteilchenfluss
    psi_e = chi_plus * exp(1i*ke_plus*(x-L)) ...
        + r_N_werte(i) .* chi_minus .* exp(1i*ke_minus*(x-L));
    psi_h = r_A_werte(i) .* chi_plus .* exp(1i*kh_plus*(x-L));
    Psi_E = [psi_e; psi_h];

    % Fluss nach Fahrplan berechnen mithilfe von Matrix A
    for ix = 1:length(x)
        Psi = Psi_E(:,ix);
        j_qp(i,ix) = real(Psi' * A_BdG * Psi);
    end

    % Konstruktion von rho_E für ein konkretes Ortspaar
    i1 = 100;
    i2 = 200;
    rho_E_werte(:,:,i) = Psi_E(:,i1) * Psi_E(:,i2)';
end


%% Ausgabe
% fprintf('rN numerisch: %.6e %+.6ei\n', real(r_N_werte), imag(r_N_werte));
% fprintf('rA numerisch: %.6f %+.6fi\n', real(r_A_werte), imag(r_A_werte));
% fprintf('rA analytisch: %.6f %+.6fi\n', real(r_A_ana), imag(r_A_ana));

% Wahrscheinlichkeiten prüfen
R_N = abs(r_N_werte).^2;
R_A = abs(r_A_werte).^2;
% fprintf('|r_N|^2 = %.6f\n', R_N);
% fprintf('|r_A|^2 = %.6f\n', R_A);
% fprintf('Summe R_N + R_A = %.12f\n', R_N+R_A);
figure
plot(E_werte/Delta0, R_A, 'LineWidth', 1.5)
hold on
plot(E_werte/Delta0, R_N, 'LineWidth', 1.5)
xlabel('E / \Delta_0')
ylabel('Reflexionswahrscheinlichkeit')
legend('R_A', 'R_N')
grid on

% Phase
figure
plot(E_werte/Delta0, angle(r_A_werte), 'LineWidth', 1.5)
hold on
% plot(E_werte/Delta0, -acos(E_werte/Delta0)-Phi_S, '--')
xlabel('E / \Delta_0')
ylabel('Phase von r_A (rad)')
legend('Numerisch')
grid on

% Quasiteilchenfluss
figure
plot(x*1e9, j_qp(100,:)/vF, 'LineWidth', 1.5)
xlabel('x (nm)')
ylabel('j_{qp} / v_F')
grid on
fprintf('Maximaler Quasiteilchenfluss: %.3e\n', max(abs(j_qp(:)/vF)));