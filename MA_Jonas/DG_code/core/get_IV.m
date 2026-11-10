function [I,J] = get_IV(p)

[SysM,rhs] = get_Diff(p);
I = zeros(p.NU,p.N_chi*p.N_K_chi);
for sz = 1:length(p.U_IV)
    fprintf('Iteration %d: %f V\n', sz, p.U_IV(sz))
    if (p.doIV == true && p.doPoisson == false)
        A = get_Drift(p,p.U_IV(sz));
        A = SysM+A;
        rho = A\rhs;
        if (p.doFEM == true)
            rho_u = (reshape(rho, p.N_chi*p.N_K_chi, p.N_xi*(p.N_K_xi-1)+1)*p.phi.');
        else
            rho_u = (reshape(rho, p.N_chi*p.N_K_chi, p.N_xi*p.N_K_xi)*p.phi.');
        end
    else
        [~,~,rho_u] = Poisson(p,p.U_IV(sz),SysM,rhs);
    end

    [~,J] = getObservables(p,rho_u);
    I(sz,:) = J;
end
    J= mean(I,2);
%     figure('name', 'IV-Kurve')
%     mesh(p.chi_Koord*1e9,U_IV,I, FaceColor='interp')
%     xlabel('\chi in nm')
%     ylabel('U in V')
%     zlabel('$J_U(x)$ in $m^{-2}$')
%     set(gca,'FontSize',20);
% 
%     figure('name', 'JV-Kurve')
    plot(p.U_IV, J)
%     xlabel('U in V')
%     ylabel('$\overline{J}_U$ in $m^{-2}$')
%     set(gca,'FontSize',20);
end
