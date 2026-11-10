function [E, V] = solve_bandstructure(valley_index, x_index, mat)

n_of_modes = mat.n_of_modes;

Ny = mat.Ny;

try mat.transv_np;
    if mat.transv_np ==1
        Ht = eval_transversal_Hamiltonian_NPB(valley_index, x_index, mat);
    elseif mat.transv_np ==2
        Ht = eval_transversal_Hamiltonian_NPB_MV(valley_index, x_index, mat);
    else
        Ht = eval_transversal_Hamiltonian(valley_index, x_index, mat);
    end
catch
    disp('catch!');
    Ht = eval_transversal_Hamiltonian(valley_index, x_index, mat);
end

[V, E] = eigs(Ht, n_of_modes, 'sm');        %smallest abs := smallest magnitude eigenvallues
E = diag(E);

[E, idx] = sort(E, 'ascend');

V = real(V(:, idx));

for IM=1:n_of_modes
    V(:,IM) = sign(real(V(end,IM)))*real(real(V(:,IM))/sqrt(sum(real(V(:,IM)').^2)));
end

E = E';
V = V';

% if x_index==1
%     figure(7);
%     hold off
%     plot(V(1,:)); hold on;
%     plot(V(2,:));
%     drawnow;
% end


% figure(),hold on;
% for IM = 1:3
%     plot(E(IM)+abs(V(IM,:)),'b');
%     plot(E(IM)*ones(13,1),'b');
% end