function [E, V] = solve_subbands_3D(mat)

%% Ändern

n_of_modes = mat.n_of_modes;
n_of_valleys = mat.n_of_valleys;

Nx = mat.Nx;
Ny = mat.Ny;
Nz = mat.Nz;

E = zeros(n_of_valleys, Nx, n_of_modes);
V = zeros(n_of_valleys, Nx, n_of_modes, Ny, Nz);

for IX = 1 : Nx
    for IV = 1 : n_of_valleys
        
        [Etemp, Vtemp] = solve_bandstructure_3D(IV, IX, mat);
        E(IV, IX, :) = Etemp;
        V(IV, IX, :, :,:) = Vtemp;
        
    end
end

%Clean up eigenvectors and eigenvalues
cleanup = true;

if cleanup == true
    Nx = mat.Nx;
    Nm = mat.n_of_modes;

    cutoff = 0.5;
    
    % FLIP DEGENERATE EIGENVECTORS FOR ALL EIGENVALUES
    for IMF = 1:Nm
        V0  = squeeze(V(1,1,IMF,:,:));
        for IX = 2:Nx
            if sum(V0.* squeeze(V(1,IX,IMF,:,:))) < -cutoff
                V(1,IX,IMF,:,:)=-V(1,IX,IMF,:,:);
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