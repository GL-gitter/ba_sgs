function G_glob = get_Drift_accel_trial_permutedShell(mat, p,U_IV, EM)

    if (p.doIV == true && p.doPoisson == false)
        p.U = U_IV;
        p.U1 = zeros(size(p.xg));
        p.U1(p.xg < p.L_chi/2+p.L_d/2 & p.xg > p.L_chi/2-p.L_d/2) = linspace(0, U_IV, sum(p.xg < p.L_chi/2+p.L_d/2 & p.xg > p.L_chi/2-p.L_d/2));
        p.U1(p.xg >= p.L_chi/2+p.L_d/2) = U_IV;
    end
    if (p.doPoisson == true && p.doIV == true)
        p.U = U_IV;
    end
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
    B = V_transient(x_drift+0.5*y_drift, mat.x*1e-9, EM, p.U, p.L_chi) ...               %% Versuch: Vorzeichenwechsel V_transient
       - V_transient(x_drift-0.5*y_drift, mat.x*1e-9, EM, p.U, p.L_chi) ...
       - 1i*get_CAP(y_drift, p.L_xi, p.delta_xi);

    if (p.doDG == true)
        W_xi = p.V_tilde_xi/p.V_xi;
        W_xiCell = repmat({W_xi},1,p.N_xi);
        W_xi = blkdiag(W_xiCell{:});
        w_xi = repmat(p.w,p.N_xi,1);
        G = zeros(p.Nqx*p.N_chi, p.N_xi*p.N_K_xi, p.N_xi*p.N_K_xi);
        for kk = 1 : p.Nqx*p.N_chi
            c_xi = B(kk,:)'.*w_xi;
            W_xi_tilde = W_xi.*c_xi;
            C = p.M_xi_glob*W_xi_tilde.'*W_xi*p.delta_xi;
            G(kk,:,:) = 1/1i*p.phi'*C*p.phi;
        end
    elseif (p.doFV == true)
        
        accel = true;
        if (accel == false)
            G = zeros(p.Nqx*p.N_chi, p.Nk, p.Nk);
            for i = 1 : p.Nqx*p.N_chi
                C = diag(B(i,1:end-1) + B(i,2:end))+diag(B(i,2:end-1),-1)+ diag(B(i,2:end-1),+1);
                C = p.delta_xi/4*C;
                G(i,:,:) = 1i*p.phi'*C*p.phi;
            end
            C = G;
        %%
        
        
        % Calculate Diags
        elseif (accel == true)      
            % Parameter und Dimensionen
            num_iterations = p.Nqx * p.N_chi;
            M_C = size(B, 2) - 1;
            Nk_val = p.Nk;
        
            % --- Schritt 1: Vorbereitung der Diagonalwerte (schleifenfrei) ---
            B_relevant = B(1:num_iterations, :);
            main_diag_elements_all = B_relevant(:, 1:M_C) + B_relevant(:, 2:M_C+1); % num_iter x M_C
            off_diag_elements_all  = B_relevant(:, 2:M_C);                   % num_iter x (M_C-1)
        
            % --- Schritt 2: Schleifenfreie Erstellung von C_3D ---
            C = zeros(M_C, M_C, num_iterations, 'like', B(1));     %C_3D
            
            % 0-basierter Seiten-Offset für lineare Indizes
            page_offset_0based = reshape(0:num_iterations-1, 1, 1, num_iterations) * (M_C * M_C);
        
            % Hauptdiagonalen
            row_col_idx_main_0based = 0:M_C-1;
            linear_idx_main_single_page_0based = row_col_idx_main_0based * M_C + row_col_idx_main_0based;
            final_linear_idx_main = 1 + bsxfun(@plus, linear_idx_main_single_page_0based, page_offset_0based); % Ergibt 1 x M_C x num_iterations
            % main_diag_elements_all ist num_iter x M_C. Permutiere zu M_C x 1 x num_iter für Zuweisung
            C(final_linear_idx_main(:)) = permute(main_diag_elements_all, [2, 1]); % Werte werden zu M_C x num_iter, dann als Vektor    %C_3D
        
            if M_C > 1
                % Untere Nebendiagonalen (diag(-1))
                row_idx_sub_0based = 1:M_C-1;    % Zeilenindizes (0-basiert)
                col_idx_sub_0based = 0:M_C-2;    % Spaltenindizes (0-basiert)
                linear_idx_sub_single_page_0based = col_idx_sub_0based * M_C + row_idx_sub_0based;
                final_linear_idx_sub = 1 + bsxfun(@plus, linear_idx_sub_single_page_0based, page_offset_0based); % 1 x (M_C-1) x num_iterations
                C(final_linear_idx_sub(:)) = permute(off_diag_elements_all, [2, 1]); % Werte (M_C-1) x num_iter, dann Vektor    %C_3D
        
                % Obere Nebendiagonalen (diag(+1))
                row_idx_super_0based = 0:M_C-2;  % Zeilenindizes (0-basiert)
                col_idx_super_0based = 1:M_C-1;  % Spaltenindizes (0-basiert)
                linear_idx_super_single_page_0based = col_idx_super_0based * M_C + row_idx_super_0based;
                final_linear_idx_super = 1 + bsxfun(@plus, linear_idx_super_single_page_0based, page_offset_0based); % 1 x (M_C-1) x num_iterations
                C(final_linear_idx_super(:)) = permute(off_diag_elements_all, [2, 1]); % Werte (M_C-1) x num_iter, dann Vektor    %C_3D
            end
            
            C = (p.delta_xi/4) * C;       %C_3D
        
            % --- Schritt 3 & 4: Matrixmultiplikationen mit pagefun und Permutation (bleiben gleich) ---
            if ~isequal(size(p.phi), [M_C, Nk_val])
                error('Dimensionen von p.phi ([%d, %d]) passen nicht zu den erwarteten Dimensionen ([%d, %d]).', ...
                      size(p.phi,1), size(p.phi,2), M_C, Nk_val);
            end
            
            C = pagemtimes(C, p.phi);    %C_3D
            C = pagemtimes(p.phi', C);   
            %Temp_3D = pagefun(@mtimes, C_3D, p.phi);
            %G_intermediate_paged = pagefun(@mtimes, p.phi', Temp_3D);
            C = 1i * permute(C, [3, 1, 2]);
            
        
            
        end
        %%
    elseif (p.doFEM == true)
        W_xi = p.V_tilde_xi/p.V_xi;
        W_xiCell = repmat({W_xi},1,p.N_xi+1);
        W_xi = blkdiag(W_xiCell{:});
        w_xi = repmat(p.w,p.N_xi+1,1);
        G = zeros(p.Nqx*p.N_chi, p.N_xi*(p.N_K_xi-1)+1, p.N_xi*(p.N_K_xi-1)+1);
        % Driftwerte in einer Tri-diagonalen Matrix aufnehmen und diagonalisieren
        for ii = 1 : p.Nqx*p.N_chi
            C_xi = zeros(p.N_xi*(p.N_K_xi-1)+1,p.N_xi*(p.N_K_xi-1)+1);
            c_xi = B(ii,:)'.*w_xi;
            W_xi_tilde = W_xi.*c_xi;
            C = W_xi_tilde.'*W_xi;
            C_xi(1:p.N_K_xi,1:p.N_K_xi)=C(1:p.N_K_xi,1:p.N_K_xi);
            for aa = 2:p.N_xi-1
                C_xi(aa:aa+p.N_K_xi-1, aa:aa+p.N_K_xi-1)=[C_xi(aa,aa)+C((aa-1)*p.N_K_xi+1,(aa-1)*p.N_K_xi+1), C((aa-1)*p.N_K_xi+1,(aa-1)*p.N_K_xi+p.N_K_xi); 
                                                          C((aa-1)*p.N_K_xi+p.N_K_xi,(aa-1)*p.N_K_xi+1),      C((aa-1)*p.N_K_xi+p.N_K_xi,(aa-1)*p.N_K_xi+p.N_K_xi)];
            end
            C_FEM = p.delta_xi/4*p.M_xi_glob*C_xi;
            G(ii,:,:) = 1/1i*p.phi'*C_FEM*p.phi;
        end
    end
    
    if (p.doDG == true)
        Ny = p.N_xi*p.N_K_xi;
    elseif (p.doFEM == true)
        Ny = p.N_xi*(p.N_K_xi-1)+1;
    elseif (p.doFV == true)
        Ny = p.Nk;
    end
    
    W = p.V_tilde/p.V;

    %c = zeros(p.Nqx*p.N_chi, Ny, Ny);
    % index = (1:p.N_K_chi * p.N_K_chi)';
    % index_zl = zeros(p.dim * p.dim / p.N_chi, 1);
    % index_sp = index_zl;
    % werte_Drift = index_zl;
    % tol = 1e-16 * norm(G(:));

    % Vorkalkulation von Konstanten außerhalb der Schleifen
    N_chi = p.N_chi;
    N_K_chi = p.N_K_chi;
    Nqx = p.Nqx;
    p.dim = N_chi*N_K_chi*Ny; 

   
    % Berechnung der gewichteten Werte mit Gauß-Lobatto-Quadratur
    C(:,:,:) = C(:, :, :) .* repmat(p.w,N_chi,1);

    W_rep    = repmat(W, N_chi, 1);
    W_rep    = repmat(W_rep, Ny*Ny,1);
    C    = reshape(C, [], 1);
    %c_rep    = reshape(c, [], 1);
    %c_rep = reshape(c_rep, [], 1);

    %W_tilde = W_rep .* c_rep;
    C = W_rep .* C;

    %C_r = reshape(C.', 3, 10, []);
    C = reshape(C.', 3, 10, []);
    W_rep   = reshape(W_rep', 3, 10, []);
    W_rep   = permute(W_rep,[2,1,3]);
    C = pagemtimes(C, W_rep);
    C = reshape(permute(C, [2,3,1]), N_K_chi*N_chi*Ny*Ny,N_K_chi);
    %W_prod = reshape(W_prod, N_K_chi*N_chi, N_K_chi*Ny*Ny);

    C = repmat((p.Q_drift * p.invM), 1, N_chi*Ny*Ny)' .* C;


    if mat.dg.params.permute_shell == false
        C = C';
        % this has to be adapted for the number of N_K_chi, try to mod for auto adaption
        %Drift1 = reshape(Drift(:, 1:3:(end-2)), [], 1);
        %Drift2 = reshape(Drift(:, 2:3:(end-1)), [], 1);
        %Drift3 = reshape(Drift(:, 3:3:end), [], 1);
        C = [reshape(C(:, 1:3:(end-2)), [], 1); reshape(C(:, 2:3:(end-1)), [], 1); reshape(C(:, 3:3:end), [], 1)];
        %werte_Drift = ;
        %index_zl(:) = zl_inner(:);
        %index_sp(:) = sp_inner(:);
    
        %zl_idx = uint32(repmat((1:(N_K_chi*N_chi*Ny))', Ny, 1));
        zl_idx = repmat((1:(N_K_chi*N_chi*Ny))', Ny, 1);

        zl_idx = repmat(zl_idx, N_K_chi,1);
    
        %sp_idx = uint32(repmat((1:N_K_chi:N_K_chi*N_chi),N_K_chi,1));
        sp_idx = repmat((1:N_K_chi:N_K_chi*N_chi),N_K_chi,1);

        sp_idx = reshape(sp_idx,N_K_chi*N_chi,1);
        sp_idx = repmat(sp_idx,Ny,1);
        sp_idx = repmat(sp_idx,Ny,1);
    
        %bias = repmat((0:N_K_chi*N_chi:(Ny-1)*N_K_chi*N_chi),N_K_chi*N_chi*Ny,1);
        %bias = reshape(bias,[],1);
        sp_idx = sp_idx + reshape(repmat((0:N_K_chi*N_chi:(Ny-1)*N_K_chi*N_chi),N_K_chi*N_chi*Ny,1),[],1);
        %sp_idx = sp_idx + bias;
        sp_idx = [sp_idx; sp_idx+1; sp_idx+2];


    elseif mat.dg.params.permute_shell == true
        %% Line Index   (N_K_chi, N_chi, Ny)
        vspace=(0:N_K_chi*Ny:(N_chi-1)*N_K_chi*Ny);
        Block=ones(N_K_chi,N_K_chi);
        
        Blockx=ones(N_chi,N_K_chi,N_K_chi);
        Blockx=Blockx.*vspace.';
        
        
        XX=[1 1 1;Ny+1 Ny+1 Ny+1;2*Ny+1 2*Ny+1 2*Ny+1];
        XX=repmat(XX,1,1,N_chi);
        XX=permute(XX, [3,1,2]);
        
        Test123=Blockx+XX;
        
        
        NN=(0:1:Ny-1).';
        NN=repmat(NN,1,N_chi);
        NN=reshape(permute(NN,[2,1]),[],1);
        
        Test123=repmat(Test123,Ny,1,1);
        
        Offset=ones(N_chi*Ny,N_K_chi,N_K_chi).*NN;
        
        Line_idx = Test123+Offset;
        
        Line_idx = permute(Line_idx, [2,1,3]);
        Line_idx = reshape(Line_idx, [],N_K_chi);
        
        Line_idx = uint32(repmat(Line_idx, Ny,1));
        %Line_idx = repmat(Line_idx, Ny,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Column Index
        
        XC=[1 Ny+1 2*Ny+1;1 Ny+1 2*Ny+1;1 Ny+1 2*Ny+1];
        XC=repmat(XC,1,1,N_chi);
        XC=permute(XC, [3,1,2]);
        
        vspace_c=(0:N_K_chi*Ny:(N_chi-1)*N_K_chi*Ny);
        
        Blockx_c=ones(N_chi,N_K_chi,N_K_chi);
        Blockx_c=Blockx_c.*vspace_c.';
        
        
        Test_c=Blockx_c+XC;
        
        Test_c=repmat(Test_c,Ny,1,1);
        Column_idx = permute(Test_c, [2,1,3]);
        Column_idx = reshape(Column_idx, [],N_K_chi);
        
        Column_idx = uint32(repmat(Column_idx, Ny,1));
        %Column_idx = repmat(Column_idx, Ny,1);
        
        spacer = uint32(repmat(ones(N_K_chi,N_K_chi),N_chi*Ny*Ny,1));
        %spacer = repmat(ones(N_K_chi,N_K_chi),N_chi*Ny*Ny,1);
        
        spacer_offset = (0:1:(Ny-1)).';
        
        %spacer_offset = repmat(spacer_offset,1,N_K_chi*N_chi*Ny);
        spacer_offset = uint32(repmat(spacer_offset,1,N_K_chi*N_chi*Ny));

        spacer_offset=reshape(permute(spacer_offset,[2,1]),[],1);
        
        spacer = spacer.*spacer_offset;
        
        Column_idx = Column_idx+spacer;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % spacer_offset = repmat(spacer_offset,1,3);
        % spacer_offset=reshape(permute(spacer_offset,[2,1]),[],1);
        
        zl_idx = reshape(Line_idx,[],1);
        
        sp_idx = reshape(Column_idx,[],1);

        C=reshape(C,[],1);
    end

    

    G_glob = [zl_idx, sp_idx, C];
    %G_glob = myspconvert(G_glob, p.dim, p.dim, 1e-15);
    G_glob = sparse(zl_idx, sp_idx, C, p.dim, p.dim);
    
end
