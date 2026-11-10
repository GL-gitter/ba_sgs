function G_j = driftmatrix_loc(p,j,werte_Drift,index_zl,index_sp,index,W,tol,Ny,G)

zl = ((1-1)*(p.N_K_chi*p.N_chi)+1 : 1*(p.N_K_chi*p.N_chi))';
    for m = 1 : Ny
        sp = (m-1)*(p.N_K_chi*p.N_chi)+1 : m*(p.N_K_chi*p.N_chi);
        for d = 1 : p.N_chi
            c(1:p.Nqx, 1) = G((d-1)*p.Nqx+1 : d*p.Nqx, j, m).*p.w;      % Gauß-Lobatto Quadratur
            if (norm(c) > tol)          % nummerische Nullen ausgrenzen
                xBlock = (d-1)*p.N_K_chi+1 : d*p.N_K_chi;
                zl_inner = zl(xBlock)*ones(1, p.N_K_chi);
                sp_inner = ones(p.N_K_chi, 1)*(sp(xBlock));
                W_tilde = W.*c;
                Drift = p.Q_drift*p.invM.*(W_tilde.'*W);
                werte_Drift(index) = Drift(:);
                index_zl(index) = zl_inner;
                index_sp(index) = sp_inner;
                index = index + p.N_K_chi*p.N_K_chi;
            end
        end
    end
G_j = [index_zl(:), index_sp(:), werte_Drift(:)];
G_j = myspconvert(G_j, p.N_chi*p.N_K_chi, p.dim, 1e-15);
end