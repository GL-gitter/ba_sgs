function G_glob = get_DriftFEM(Q_drift, delta_chi, delta_xi, L_xi, Nqx, N_chi, N_Knoten, N_xi, N_K_xi, disk_Rechengebiet1, phi, t_local, U0, U, rampTime, V0, dim, M_xi, M, invM, xi_Rechengebiet)
    
    M_xi_glob = zeros(N_xi*(N_K_xi-1)+1,N_xi*(N_K_xi-1)+1);
    M_xi_glob(1:N_K_xi,1:N_K_xi)=M_xi;
    for aa = 1:N_xi-1
        M_xi_glob(aa+1:aa+N_K_xi, aa+1:aa+N_K_xi)=[M_xi_glob(aa+1,aa+1)+M_xi(1,1), M_xi(1,N_K_xi); 
                                                    M_xi(N_K_xi,1), M_xi(N_K_xi,N_K_xi)];
    end
%     M_xi_glob = M_xi_glob(2:end-1, 2:end-1);

    % Gauß-Lobatto Kantenkoordinaten bestimmen
    [x, w] = JacobiGL(0, 0, Nqx-1);
    r = JacobiGL(0, 0, N_Knoten-1);
    r_xi = JacobiGL(0, 0, N_K_xi-1);
    V_tilde = Vandermonde1D(N_Knoten-1, x);
    V = Vandermonde1D(N_Knoten-1, r);
    V_tilde_xi = Vandermonde1D(N_K_xi-1, x);
    V_xi = Vandermonde1D(N_K_xi-1, r_xi);
    x_anteil = diff(x);                 % Abstände zwischen den GL-Punkten
    x_prozent = x_anteil/2;             % Prozentualer Anteil der Quadraturpunkte
    chi_anteil = x_prozent * delta_chi;
    xi_anteil = x_prozent * delta_xi;
    chi_GL_koord = zeros(1, Nqx);
    xi_GL_koord = zeros(1,Nqx);
    
    % Quadratur-Koordinaten für ein Chi-Element bestimmen
    for i = 1 : Nqx-1
        chi_GL_koord(i+1) = chi_GL_koord(i)+chi_anteil(i);
    end

    % Quadratur-Koordinaten für ein Xi-Element bestimmen
    for i = 1:Nqx-1
        xi_GL_koord(i+1) = xi_GL_koord(i)+xi_anteil(i);
    end
    
    % Quadratur-Koordinaten für alle Chi-Elemente bestimmen
    Kantenkoord_GL = zeros(1, N_chi*Nqx);
    for i = 0 : N_chi-1
        Kantenkoord_GL(1, i*Nqx+1 : (i+1)*Nqx) = disk_Rechengebiet1(i+1) + chi_GL_koord;
    end

    % Quadratur-Koordinaten für alle Xi-Elemente bestimmen
    xiKantenkoord_GL = zeros(1, N_xi*Nqx);
    for i = 0 : N_xi
        xiKantenkoord_GL(1, i*Nqx+1:(i+1)*Nqx) = xi_Rechengebiet(i+1)+xi_GL_koord;
    end
    
    % chi-Koordinaten für den Driftoperator
    x_drift = zeros(N_chi*Nqx, N_xi*(N_K_xi-1)+1);
    for i = 1 : (N_xi+1)*Nqx
        x_drift(:, i) = Kantenkoord_GL.';
    end
    x_drift_sym = x_drift-max(x_drift)/2;
    
    % xi-Koordinaten für den Driftoperator
    y_drift = zeros(N_chi*Nqx, (N_xi+1)*Nqx);
    for i = 1 : Nqx*N_chi
        y_drift(i, :) = xiKantenkoord_GL;
    end
    
    % k = -2*pi/L_xi*((1:N_xi)-1/2*(N_xi+1));
    % dk = abs(k(2)-k(1));
    
    % Driftwerte bestimmen und in einer Matrix (Nqx*N_chi, N_xi, N_xi) aufnehmen
    B = Potential(x_drift+0.5*y_drift, N_xi*(N_K_xi-1)+1) - V_transient_old(x_drift_sym+0.5*y_drift, t_local, U0, U, rampTime) ...
      - Potential(x_drift-0.5*y_drift, N_xi*(N_K_xi-1)+1) + V_transient_old(x_drift_sym-0.5*y_drift, t_local, U0, U, rampTime) ...
      - 1i*get_CAP(y_drift, L_xi, delta_xi)/V0;
    
    G = zeros(Nqx*N_chi, N_xi*(N_K_xi-1)+1, N_xi*(N_K_xi-1)+1);
    
    W_xi = V_tilde_xi/V_xi;
    W_xiCell = repmat({W_xi},1,N_xi+1);
    W_xi = blkdiag(W_xiCell{:});
    w_xi = repmat(w,N_xi+1,1);

    % Driftwerte in einer Tri-diagonalen Matrix aufnehmen und diagonalisieren
    for ii = 1 : Nqx*N_chi
        C_xi = zeros(N_xi*(N_K_xi-1)+1,N_xi*(N_K_xi-1)+1);
        c_xi = B(ii,:)'.*w_xi;
        W_xi_tilde = W_xi.*c_xi;
        C = W_xi_tilde.'*W_xi;
        C_xi(1:N_K_xi,1:N_K_xi)=C(1:N_K_xi,1:N_K_xi);
        for aa = 2:N_xi-1
            C_xi(aa:aa+N_K_xi-1, aa:aa+N_K_xi-1)=[C_xi(aa,aa)+C((aa-1)*N_K_xi+1,(aa-1)*N_K_xi+1), C((aa-1)*N_K_xi+1,(aa-1)*N_K_xi+N_K_xi); 
                                                  C((aa-1)*N_K_xi+N_K_xi,(aa-1)*N_K_xi+1), C((aa-1)*N_K_xi+N_K_xi,(aa-1)*N_K_xi+N_K_xi)];
        end
        C_FEM = C_xi*M_xi_glob;
        G(ii,:,:) = 1/1i*phi'*C_FEM*phi;
    end
    
    W = V_tilde/V;
    c = zeros(Nqx, 1);
    index = (1:N_Knoten*N_Knoten)';
    index_zl = zeros(dim*dim/N_chi,1);
    index_sp = index_zl;
    werte_Drift = index_zl;
    tol = 1e-16*norm(G(:));
    
    for j = 1 : N_xi*(N_K_xi-1)
        zl = ((j-1)*(N_Knoten*N_chi)+1 : j*(N_Knoten*N_chi))';
        for m = 1 : N_xi*(N_K_xi-1)
            sp = (m-1)*(N_Knoten*N_chi)+1 : m*(N_Knoten*N_chi);
            for d = 1 : N_chi
                c(1:Nqx, 1) = G((d-1)*Nqx+1 : d*Nqx, j, m).*w;      % Gauß-Lobatto Quadratur       
                
                if (norm(c) > tol)          % nummerische Nullen ausgrenzen
                    xBlock = (d-1)*N_Knoten+1 : d*N_Knoten;
                    zl_inner = zl(xBlock)*ones(1, N_Knoten);
                    sp_inner = ones(N_Knoten, 1)*(sp(xBlock));
                    W_tilde = W.*c;
                    Drift = Q_drift*delta_xi/12*invM*W_tilde.'*W;
                    werte_Drift(index) = Drift(:);
                    index_zl(index) = zl_inner;
                    index_sp(index) = sp_inner;
                    index = index + N_Knoten*N_Knoten;
                end
            end
        end
    end
    
    G_glob = [index_zl(:), index_sp(:), werte_Drift(:)];
    G_glob = myspconvert(G_glob, dim, dim, 1e-15);       % Aufstellen der glob. Driftmatrix und Ignorieren nummerischer Nullen (1e-15)

end