function [A,rhs,rho_u] = Poisson(p,U_IV,GDI,rhs)

    k = 0;
    if (p.doPoisson == true && p.doIV == false)
        [GDI, rhs] = get_Diff(p);
        U_IV = p.U;
    end 

    if (p.doPoisson == true && p.doIV == true)
        p.U1 = zeros(size(p.xg));
        p.U1(p.xg < p.L_chi/2+p.L_d/2 & p.xg > p.L_chi/2-p.L_d/2) = linspace(0, U_IV, sum(p.xg < p.L_chi/2+p.L_d/2 & p.xg > p.L_chi/2-p.L_d/2));
        p.U1(p.xg >= p.L_chi/2+p.L_d/2) = U_IV;
        p.UU = p.U1;
    end

    while p.err > 1e-3
        p.U1 = 0.7*p.U1+0.3*p.UU;

        A = get_Drift(p,U_IV);

        A = GDI+A;
        rho = A\rhs;
        if (p.doFEM == true)
            rho_u = (reshape(rho, p.N_chi*p.N_K_chi, p.N_xi*(p.N_K_xi-1)+1)*p.phi.');
        else
            rho_u = (reshape(rho, p.N_chi*p.N_K_chi, p.N_xi*p.N_K_xi)*p.phi.');
        end
        [~,n] = getN(p,rho_u,p.N_K_chi+1);
        
        P = diag(-(p.eps_V(1:p.Np)+p.eps_V(2:p.Np+1))./(p.dx.^2)-p.q*n./p.V_ref, 0) + ...
            diag(  p.eps_V(2:p.Np)./(p.dx.^2),   -1) + ...
            diag(  p.eps_V(2:p.Np)./(p.dx.^2),   +1);
        R = p.q*(n.*(1-p.U1'./p.V_ref)-p.Nd_V).';
        R(p.Np) = R(p.Np) - U_IV*P(p.Np-1, p.Np);
        p.UU = (P\R);

        p.err = max(max(abs(p.UU-p.U1)));
        k = k+1;
        fprintf('Fehler %d-te Iteration: %f\n',k,p.err);
 
    end

%% Plot
if (p.doIV == false)
    if (p.doFV == true || p.doFEM == true)
        figure('name','Self-consistent real part')
        mesh(p.chi_Koord, p.xi_Rechengebiet,real(rho_u'), FaceColor="interp")
    
        figure('name','Self-consistent imaginary part')
        mesh(p.chi_Koord, p.xi_Rechengebiet,imag(rho_u'), FaceColor="interp")
    else
        figure('name','Self-consistent real part')
        mesh(p.chi_Koord,p.xi_Koord,real(rho_u'), FaceColor='interp')

        figure('name','Self-consistent imaginary part')
        mesh(p.chi_Koord,p.xi_Koord,imag(rho_u'), FaceColor='interp')
    end
end
end
