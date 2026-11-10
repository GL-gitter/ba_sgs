function [n,J] = getObservables(p,rho)

%     if (p.doDG == true || p.doFV == true)
%         Dy = spdiags(ones(p.N_xi*p.N_K_xi,1)*[-1/2/p.delta_xi, 1/2/p.delta_xi], [-1,1], p.N_xi*p.N_K_xi,p.N_xi*p.N_K_xi);
%     elseif (p.doFEM == true)
%         Dy = spdiags(ones((p.N_K_xi-1)*p.N_xi+1,1)*[-1/2/p.delta_xi, 1/2/p.delta_xi], [-1,1], (p.N_K_xi-1)*p.N_xi+1,(p.N_K_xi-1)*p.N_xi+1);
%     end
    
    rho_imag = imag(rho)*p.Dy;
    J0 = p.hquer/p.m_elektron*p.q;

    if (p.doDG == true)
        J = J0*(rho_imag(:,p.N_K_xi*p.N_xi/2-2)+rho_imag(:,p.N_K_xi*p.N_xi/2+2))/2;
        n = (real(rho(:,p.N_K_xi*p.N_xi/2+1))+real(rho(:,p.N_K_xi*p.N_xi/2)))/2;
    elseif ( p.doFV == true)
        J = -J0*(rho_imag(:,p.N_xi/2+1)+rho_imag(:,p.N_xi/2))/2;
        n = (real(rho(:,p.N_xi/2+1))+real(rho(:,p.N_xi/2)))/2;
    elseif (p.doFEM == true)
        J = J0*(rho_imag(:,((p.N_K_xi-1)*p.N_xi+1)/2+2)+rho_imag(:,((p.N_K_xi-1)*p.N_xi+1)/2-2))/2;
        n = (real(rho(:,((p.N_K_xi-1)*p.N_xi+1)/2+1))+real(rho(:,((p.N_K_xi-1)*p.N_xi+1)/2)))/2;
    end
end