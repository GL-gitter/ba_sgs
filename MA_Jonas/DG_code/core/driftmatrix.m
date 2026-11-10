function G_j = driftmatrix(p,j,Ny,W_all,G,w_all,invM_all,tol)
    G_j = zeros(p.N_K_chi*p.N_chi,p.dim);
    for m = 1:Ny
        c = G(:,j,m).*w_all;
        if (norm(c) > tol)
            W_tilde = W_all.*c;
            G_j(:,(m-1)*p.N_chi*p.N_K_chi+1:m*p.N_chi*p.N_K_chi) = p.Q_drift*invM_all.*(W_tilde.'*W_all);
        end
    end
end