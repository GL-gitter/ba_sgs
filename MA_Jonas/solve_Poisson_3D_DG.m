function [Vbi, iter_err] = solve_Poisson_3D_DG(M, node_type, node_num, n, mat, Vbi, Vg_bias)

eps0 = 8.854E-12;
q   = 1.602E-19;
kB  = 1.38E-23;
T   = mat.Temp;
N = mat.Nd*1E+6;
alpha = mat.poisson.opt.alpha;

Vbi_old = Vbi;

Nx = mat.Nx;
Ny = mat.Ny;
Nz = mat.Nz;

%% Adaption for DG algorithm
N_K_chi = mat.dg.params.N_K_chi;
N_chi = ((mat.Nx-1)/(N_K_chi-1));
Np = N_chi*N_K_chi;
%%
node_cnt = 0;

free_node = node_num(node_type==1);
bnd_node = node_num(node_type==2);


for IX = 1 : Np
    for IY = 1 : Ny
        for IZ = 1 : Nz
        node_cnt = node_cnt+1;
        
            n_p(node_cnt,1) = n(IX, IY, IZ);
            g(node_cnt,1)   = Vg_bias*mat.boundary(IX, IY, IZ);
            N_p(node_cnt,1) = N(IX, IY, IZ);
            
            Vbi_old_p(node_cnt,1) = Vbi_old(IX, IY, IZ);
        end
    end
    
end



switch mat.poisson.opt.solve
    
    case 'direct'
        Vbi_new = zeros(node_cnt,1);
        Vbi_new(bnd_node) = g(bnd_node);
        RHS = -q*(N_p-n_p)+M*Vbi_new;
        
        Vbi_new(free_node) = -M(free_node,free_node)\RHS(free_node);
        Vbi_new = Vbi_old_p + alpha*(Vbi_new-Vbi_old_p);
        iter_err = norm(Vbi_new-Vbi_old_p, 2);
        node_cnt = 0;
        
        for IX = 1 : Np
            for IY = 1 : Ny
                for IZ = 1 : Nz
                    node_cnt = node_cnt+1;
                    Vbi(IX,IY,IZ) = +Vbi_new(node_cnt);
                end
            end
        end
        
    case 'newton-rhapson'
        
        Vbi_new = zeros(node_cnt,1);
        Vbi_old_p(bnd_node) = g(bnd_node);
        Vbi_new(bnd_node) = g(bnd_node);
        
        
        dVbi = zeros(node_cnt,1);
        
        F     = M*Vbi_old_p - q*(N_p-n_p);
                 % J     = (M-q^2*spdiags(n_p./(kB*T)));
        J     = M-spdiags(q^2*n_p./(kB*T),0,Np*Ny*Nz,Np*Ny*Nz);

        dVbi(free_node)  = (-J(free_node,free_node)\F(free_node)); % first minus due to -q*PHI = V
        Vbi_new(free_node) = Vbi_old_p(free_node)+alpha*dVbi(free_node);
        iter_err = norm(dVbi, 2);%/node_cnt;
        node_cnt = 0;
        for IX = 1 : Np
            for IY = 1 : Ny
                for IZ = 1 : Nz
                    node_cnt = node_cnt+1;
                    Vbi(IX,IY,IZ) = Vbi_new(node_cnt);
                end
            end
        end
end