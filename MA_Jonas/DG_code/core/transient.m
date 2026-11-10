function [n_over_t, j_over_t, rho_T] = transient(p, rhs, A, rho, rho_U)

n_over_t = zeros(p.N_chi*p.N_K_chi, p.tsteps);
j_over_t = n_over_t;
fid = figure('name', "Dichte");
%% Transiente Berechnung mit Runge-Kutta 4. Ordnung
for tstep = 1 : p.tsteps
    if mod((tstep-1),50) == 0
        display(['Zeitschritt ', num2str(tstep), '/', num2str(p.tsteps)])
    end
    % Runge-Kutta 4. Ordnung
    k1 = rhs-A*rho;
    k2 = rhs-A*(rho+p.dt/2*k1);
    k3 = rhs-A*(rho+p.dt/2*k2);
    k4 = rhs-A*(rho+p.dt*k3);
    rho = rho+p.dt/6*(k1+2*k2+2*k3+k4);

    result_matrix = reshape(rho, p.N_chi*p.N_K_chi, p.N_xi*p.N_K_xi)*p.phi.';
    
    % Observablen speichern
    n_over_t(:, tstep) = real(result_matrix(:, p.N_xi*p.N_K_xi/2+1));   
    j_over_t(:, tstep) = p.q*p.hquer./p.m_elektron.*imag(result_matrix(:, p.N_xi*p.N_K_xi/2+3)-result_matrix(:, p.N_xi*p.N_K_xi/2))/2/p.delta_xi; 
    % Dichteverlauf plotten
    if mod((tstep-1),5) == 0
        figure(fid)
        title(['t=', num2str(p.t), 'fs']);
        plot(p.chi_Koord,real(rho_U(p.N_xi*p.N_K_xi/2,:)), 'b');
        hold on 
        plot(p.chi_Koord,real(result_matrix(:, p.N_xi*p.N_K_xi/2)), 'r');
        hold off
    end
end

n_over_t = n_over_t(:, 1:10:end);
j_over_t = j_over_t(:, 1:10:end);
rho_T = reshape(rho, p.N_chi*p.N_K_chi, p.N_xi*p.N_K_xi)*p.phi.';