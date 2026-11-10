function [E, Vout] = solve_bandstructure_3D(valley_index, x_index, mat)

n_of_modes = mat.n_of_modes;

Ny = mat.Ny;
Nz = mat.Nz;

Ht = eval_transversal_Hamiltonian_3D(valley_index, x_index, mat);


[V, E] = eigs(Ht, n_of_modes, 'sm');
E = diag(E);

[E, idx] = sort(E, 'ascend');
E = E';
V = real(V(:, idx));
Vout = zeros(n_of_modes,Ny,Nz);
for IM=1:n_of_modes
    Vtemp = real(real(V(:,IM))/sqrt(sum(real(V(:,IM)').^2)));
    Vtemp = sign(Vtemp(1))*reshape(Vtemp,[ Ny, Nz]);

    Vout(IM,:,:)=Vtemp;
end




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