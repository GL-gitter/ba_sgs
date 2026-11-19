% Erst "Schritt3.m" ausführen für Variable sigma_x, phi_werte, E_ABS_plus_num,
% E_ABS_minus_num, k_e, k_h, M_L im Workspace

%% Definitionen
E_ABS = [E_ABS_plus_num(:), E_ABS_minus_num(:)];
x_werte = linspace(0,L,100);
Fehler_1 = 0;
Fehler_2 = 0;
Fehler_3 = 0;
Fehler_4 = 0;
Fehler_C = 0;

%% Berechnung der Eigenfunktionen
psi_E_n_x = cell(2,length(phi_werte));
for j = 1:length(phi_werte)
    phi = phi_werte(j);
    for n = 1:2
        E_n = E_ABS(j,n);
        C = T_h(E_n)*M_L(E_n,phi) - M_R(E_n,phi)*T_e(E_n);
        % b = -C(1,1)/C(1,2);
        % psi_e_0 = [1; b];
        [~,~,V] = svd(C);
        psi_e_0 = V(:,end);
        Fehler_C = max(Fehler_C, norm(C*psi_e_0));
        psi_e_x = @(x) (cos(k_e(E_n)*x)*eye(2)+1i*sin(k_e(E_n)*x)* ...
                        sigma_x)*psi_e_0;
        psi_h_x = @(x) (cos(k_h(E_n)*x)*eye(2)+1i*sin(k_h(E_n)*x)*sigma_x)* ...
                        M_L(E_n,phi)*psi_e_0;
        psi_E_n_x{n,j} = @(x) [psi_e_x(x); psi_h_x(x)];

        % Normierung über das Graphengebiet
        psi = psi_E_n_x{n,j};

        N = sqrt(integral(@(x) sum(abs(psi(x)).^2), 0, L, 'ArrayValued', true));

        psi_E_n_x_normed{n,j} = @(x) psi(x)/N;
    end
end
fprintf('Maximaler Fehler der Eigenwertbedingung: %.3e\n', Fehler_C);

%% Dichtematrix rho_E konstruieren
rho_E = cell(2,length(phi_werte));
for j = 1:length(phi_werte)
    for n = 1:2
        psi = psi_E_n_x_normed{n,j};
        rho_E{n,j} = @(x1,x2) psi(x1)*psi(x2)';
    end
end

%% Überprüfen der Andreev-RB
for j = 1:length(phi_werte)
    phi = phi_werte(j);
    for n = 1:2
        E_n = E_ABS(j,n);
        ML = M_L(E_n,phi);
        MR = M_R(E_n,phi);
        for ix = 1:length(x_werte)
            % Werte aus der Dichtematrix holen
            x = x_werte(ix);
            rho_0_x = rho_E{n,j}(0,x);
            rho_x_0 = rho_E{n,j}(x,0);
            rho_L_x = rho_E{n,j}(L,x);
            rho_x_L = rho_E{n,j}(x,L);
            A1 = rho_0_x(3:4,1:2);
            B1 = ML*rho_0_x(1:2,1:2);
            A2 = rho_x_0(1:2,3:4);
            B2 = rho_x_0(1:2,1:2)*ML';
            A3 = rho_L_x(3:4,1:2);
            B3 = MR*rho_L_x(1:2,1:2);
            A4 = rho_x_L(1:2,3:4);
            B4 = rho_x_L(1:2,1:2)*MR';
            % Vergleich mit Andreev-RB über Frobeniusnorm
            Fehler_1 = max(Fehler_1, norm(A1-B1,'fro')/ ...
                           (norm(A1,'fro')+norm(B1,'fro')));
            Fehler_2 = max(Fehler_2, norm(A2-B2,'fro')/ ...
                           (norm(A2,'fro')+norm(B2,'fro')));
            Fehler_3 = max(Fehler_3, norm(A3-B3,'fro')/ ...
                           (norm(A3,'fro')+norm(B3,'fro')));
            Fehler_4 = max(Fehler_4, norm(A4-B4,'fro')/ ...
                           (norm(A4,'fro')+norm(B4,'fro')));
        end
    end
end

% Ausgabe der Fehlerwerte
fprintf('Abweichung bei rho_he(0,x2): %.3e\n', Fehler_1);
fprintf('Abweichung bei rho_eh(x1,0): %.3e\n', Fehler_2);
fprintf('Abweichung bei rho_he(L,x2): %.3e\n', Fehler_3);
fprintf('Abweichung bei rho_he(x1,L): %.3e\n', Fehler_4);