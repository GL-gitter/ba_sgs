function [Ea, Va] = solve_subbands_CMS(mat)


%% UNPACK VARIABLES
n_of_modes = mat.n_of_modes;
n_of_valleys = mat.n_of_valleys;

Nx = mat.Nx;
Ny = mat.Ny;
 
Ea = zeros(n_of_valleys, Nx, n_of_modes);
Va = zeros(n_of_valleys, Nx, n_of_modes, Ny);


%% SOLVE 2D HAMILTONIAN AND EXTRACT EV
Vtempa = zeros(n_of_modes,Ny,Nx);
Etempa = zeros(n_of_modes,Nx);


wdw=25; %5;
for IX=1+floor(wdw/2):Nx-floor(wdw/2)
    if IX==80
        disp('f');
    end
    Ht = eval_transversal_Hamiltonian_temp_new_idee(1, IX-floor(wdw/2):IX+floor(wdw/2), mat);
    [V, E] = eigs(Ht, n_of_modes, 'sm');
    E=diag(E);
    [E, idx] = sort(E, 'ascend');
    Etempa(:,IX)=E;
    V = real(V(:, idx));
    for IM = 1:n_of_modes
        Vtemp = reshape(V(:,IM),[ Ny, wdw]);
        % Vtemp = reshape(V(:,IM),[ Ny, wdw*2-1]);
        % if IX>77 && IX < 120
        %     figure(76),imagesc(abs(Vtemp)-0.005*mat.V(IX-floor(wdw/2):IX+floor(wdw/2),:)');
        % end
        Vtempa(IM,:,IX) = Vtemp(:,ceil(wdw/2));
    end
end

for IX=1:floor(wdw/2)
    Vtempa(:,:,IX)          = Vtempa(:,:,floor(wdw/2)+1);
    Vtempa(:,:,Nx+1-IX)     = Vtempa(:,:,Nx-floor(wdw/2));
    
    Etempa(:,IX)            = Etempa(:,floor(wdw/2)+1);
    Etempa(:,Nx+1-IX)       = Etempa(:,Nx-floor(wdw/2));
end


for IM = 1:n_of_modes
    Ea(1, :, IM)=squeeze(sgolayfilt(Etempa(IM,:),0,3));
end


%% MAKE SURE THE EV SIGNS ARE CORRECT AND NORMALIZE EV
for IM=1:n_of_modes
    for IX=1:Nx
        Va(1,IX,IM,:) = sign(real(Vtempa(IM,end,IX)))*real(Vtempa(IM,:,IX))/sqrt(sum(real(Vtempa(IM,:,IX)).^2,'all'));
    end
end

% figure(),imagesc(squeeze(Va(1,:,3,:)));


%% SOLVE 1D HAMILTONIAN TO OBTAIN SUBBAND ENERGIES

% for IM=1:n_of_modes
%     for IX=1:Nx
%         Vtemp= squeeze(Va(1,IX,IM,:));
%         Ht = eval_transversal_Hamiltonian(1, IX, mat);
%         Ephi = (Ht*Vtemp)./Vtemp;
%         Etempa(IM,IX)=Ephi(13);
%     end
% end

%Clean up eigenvectors and eigenvalues
cleanup = false;

if cleanup == true
    Nx = mat.Nx;
    Nm = mat.n_of_modes;

    cutoff = 0.5;
    
    % FLIP DEGENERATE EIGENVECTORS FOR ALL EIGENVALUES
    for IMF = 1:Nm
        V0  = squeeze(Va(1,1,IMF,:));
        for IX = 2:Nx
            if sum(V0.* squeeze(Va(1,IX,IMF,:))) < -cutoff
                Va(1,IX,IMF,:)=-Va(1,IX,IMF,:);
            end

        end
    end
    % 
    % % GET DEGENERATE EIGENVALUES
    % nhisd = zeros(Nm,1); %next higher is degenerate
    % for IM = 1:Nm-1
    %     nhisd(IM) = abs(EM(1,1,IM)-EM(1,1,IM+1))<0.1;
    % end
    % idx = find(nhisd>0.5);
    % 
    % % SWAP DEGENERATE EIGENVECTORS
    % for IMS = idx'
    %     % Lower Mode
    %     V0  = squeeze(mat.VSM(IMS,:,:));                                %First Slice of Lower Mode
    %     V0_const = repmat(V0,Nx,1);
    %     id_help = sum(squeeze(VM(IV,:,IMS,:,:)).*V0_const, [1,2]);
    %     idx_swap_lower1 = -0.1 < id_help;
    %     idx_swap_lower2 = id_help < cutoff;
    %     idx_swap_lower=idx_swap_lower1 & idx_swap_lower2;
    % 
    %     % Upper Mode
    %     V0  = squeeze(mat.VSM(IMS+1,:,:));                              %First Slice of Upper Mode
    %     V0_const = repmat(V0,Nx,1);
    %     id_help = sum(squeeze(VM(IV,:,IMS+1,:,:)).*V0_const, [1,2]);
    %     idx_swap_upper1 = -0.1 < id_help;
    %     idx_swap_upper2 = id_help < cutoff;
    %     idx_swap_upper=idx_swap_upper1 & idx_swap_upper2;
    % 
    %     idx_swap = find(idx_swap_lower & idx_swap_upper);
    %     %disp([IMS,'Swapping',idx_swap])
    % 
    %     Vtemp = zeros(size(VM(IV,:,IMS,:,:)));
    %     Vtemp(idx_swap,:,:) = VM(IV,idx_swap,IMS,:,:);
    % 
    %     VM(IV,idx_swap,IMS,:,:)     = VM(IV,idx_swap,IMS+1,:,:);
    %     VM(IV,idx_swap,IMS+1,:,:)   = Vtemp(idx_swap,:,:);
    % 
    %     Etemp = zeros(size(EM(IV,:,IMS)));
    %     Etemp(idx_swap) = EM(IV,idx_swap,IMS);
    % 
    %     EM(IV,idx_swap,IMS)     = EM(IV,idx_swap,IMS+1);
    %     EM(IV,idx_swap,IMS+1)   = Etemp(idx_swap);
    % end
end
end