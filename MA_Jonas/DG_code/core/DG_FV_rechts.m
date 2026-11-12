%% Lokale Funktion: rechte Seite für alle DG-Elemente und FV-Zellen
function dUdt = DG_FV_rechts(U, S, M_k, d_xi, A_X, A_X_plus, A_X_minus, ...
    A_xi_plus, A_xi_minus, S_pot) % D_F, hbar, xi für lineares Potential 
                                  % (anstelle von S_pot
    K = size(U,3);
    N_xi = size(U,4);
    dUdt = zeros(size(U));
    % Schleife über alle DG-Elemente
    for k = 1:K
        % Schleife über FV-Zellen
        for j = 1:N_xi
            j_links = mod(j-2,N_xi)+1; % periodischer linker Nachbar
            j_rechts = mod(j,N_xi)+1; % periodischer rechter Nachbar
    
            U_mitte = U(:,:,k,j);  % 4x2-Matrix, die die Dichtematrixwerte 
            % eines DG-Elements zu eier FV-Zelle
            % enthält
    
            % DG-Randfluesse in X
            if k==1
                F_links = A_X_minus*U_mitte(:,1);
            else
                F_links = A_X_minus*U_mitte(:,1) + A_X_plus*U(:,2,k-1,j);
            end
            if k==K
                F_rechts = A_X_plus*U_mitte(:,2);
            else
                F_rechts = A_X_plus*U_mitte(:,2) + A_X_minus*U(:,1,k+1,j);
            end
    
            % Schwache DG-Formulierung
            F = A_X*U_mitte;
            linke_Seite = (F*S - F_rechts*[0 1] + F_links*[1 0])/M_k;
    
            % FV-Flussdifferenz in xi
            G = (A_xi_plus*(U_mitte-U(:,:,k,j_links)) ...
                + A_xi_minus*(U(:,:,k,j_rechts)-U_mitte))/d_xi;
    
            % Potentialterm
            S_pot_term = U_mitte .* reshape(S_pot(:,k,j),1,2);
            % S_pot = D_F*xi(j)/(1i*hbar) * U_mitte;
    
            % Zeitliche Ableitung: dU/dt = DG - G
            dUdt(:,:,k,j) = linke_Seite-G + S_pot_term;
        end
    end
end