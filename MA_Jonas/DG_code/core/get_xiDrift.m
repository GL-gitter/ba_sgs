function G = get_xiDrift(p)
    p.M = p.delta_xi/2*p.M;
    p.invM = inv(p.M);
    
    %% Gauß-Lobatto coordinates
    x_anteil = diff(p.x);                 % Abstände zwischen den GL-Punkten
    x_prozent = x_anteil/2;             % Prozentualer Anteil der Quadraturpunkte
    chi_anteil = x_prozent * p.delta_chi;
    chi_GL_koord = zeros(1, p.Nqx);
    xi_GL_koord = zeros(1, p.Nqx);
    xi_anteil = x_prozent * p.delta_xi;
    % Quadratur-Koordinaten für ein Xi-Element bestimmen
    for i = 1:p.Nqx-1
        xi_GL_koord(i+1) = xi_GL_koord(i)+xi_anteil(i);
    end
    % Quadratur-Koordinaten für alle Xi-Elemente bestimmen
    if (p.doFV == true)
        xiKantenkoord_GL = linspace(-p.L_xi/2, p.L_xi/2, p.N_xi+1);
    elseif (p.doDG == true)
        xiKantenkoord_GL = zeros(1, p.Ny);
        for i = 0 : p.N_xi-1
            xiKantenkoord_GL(1, i*p.Nqx+1:(i+1)*p.Nqx) = p.xi_Rechengebiet(i+1)+xi_GL_koord;
        end
    elseif (p.doFEM == true)
        xiKantenkoord_GL = zeros(1, p.Ny);
        for i = 0 : p.N_xi
            xiKantenkoord_GL(1, i*p.Nqx+1:(i+1)*p.Nqx) = p.xi_Rechengebiet(i+1)+xi_GL_koord;
        end
    end
    % xi-Koordinaten für den Driftoperator
    y_drift = zeros(p.N_chi*p.Nqx, p.Ny);
    for i = 1 : p.Nqx*p.N_chi
        y_drift(i, :) = xiKantenkoord_GL;
    end

    % Quadratur-Koordinaten für ein Chi-Element bestimmen
    for i = 1:p.Nqx-1
        chi_GL_koord(i+1) = chi_GL_koord(i)+chi_anteil(i);
    end
    
    % Quadratur-Koordinaten für alle Chi-Elemente bestimmen
    Kantenkoord_GL = zeros(1, p.N_chi*p.Nqx);
    for i = 0 : p.N_chi-1
        Kantenkoord_GL(1, i*p.Nqx+1 : (i+1)*p.Nqx) = p.disk_Rechengebiet1(i+1) + chi_GL_koord;
    end

    % chi-Koordinaten für den Driftoperator
    x_drift = zeros(p.N_chi*p.Nqx, p.Ny);
    for i = 1 : p.Ny
        x_drift(:, i) = Kantenkoord_GL.';
    end
    
    % Driftwerte bestimmen und in einer Matrix (p.Nqx*p.N_chi, p.N_xi, p.N_xi) aufnehmen
    B = Potential(x_drift+0.5*y_drift, p) - V_transient(x_drift+0.5*y_drift, p.xg, p.U1, p.U, p.L_chi) ...
      - Potential(x_drift-0.5*y_drift, p) + V_transient(x_drift-0.5*y_drift, p.xg, p.U1, p.U, p.L_chi) ...
      - 1i*get_CAP(y_drift, p.L_xi, p.delta_xi);
    N = p.N_chi*p.N_xi;
    % G = zeros(N);
    i =1;
        C = diag(B(i,1:end-1) + B(i,2:end))+diag(B(i,2:end-1),-1)+ diag(B(i,2:end-1),+1);
        C = p.delta_xi/4*C;
        G = C;
    
    

end