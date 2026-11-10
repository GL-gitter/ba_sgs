function [Dy,R,D,k,dk,Q] = Basisfunktionen_Q_m(p, N_xi)    %[p.Dy, p.phi,p.D,p.k,p.dk,p.Q]

    if (p.doFV == true)                                 %%%% EVs and Eigenvalues for Basis transformation according to L.Schulz to implement inflow boundary conditions
        Kpy = N_xi+1;
        a = 0;
        b = +1i/2;c = -1i/2; factor = sqrt(b*c);

        eigs_A = zeros(Kpy-1,1);
        R = zeros(Kpy-1);
        norm_R=0;

        k=1:N_xi;
        for j=1:Kpy-1
            eigs_A(j) = (a - sign(imag(c)) * 2*factor*cos(pi*j/((Kpy-1)+1)));
            R(:,j) = (b/c).^(k/2).*sin(k*pi*j/Kpy);
            norm_R = norm_R+R(1,j)^2;
        end
        R = 1/sqrt(norm_R) * R;
        D = diag(eigs_A);
        Dy = spdiags(ones(p.N_xi,1)*[-1/2/p.delta_xi, 1/2/p.delta_xi], [-1,1], p.N_xi,p.N_xi);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        y = p.xi_Rechengebiet_Elem;
        Q = zeros(p.N_xi, p.Nk);
            
        k     = -2*pi/p.L_xi*((1:p.Nk)-1/2*(p.Nk+1));
        dk    = abs(k(2)-k(1));
        
        for IN = 1 : p.Nk
            Q(:,IN) = 1./sqrt(p.L_xi)*exp(+1i*k(IN)*y');
        end
        k = diag(k);
        D = k;
        R=Q;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif (p.doDG == true)
        %% Does plane wave approach can easily work here?
        if (strcmpi(p.BT,'EIG')==true)
            Dy = get_Discretization(p);
            [phi, Wellenintervall] = eigs(Dy, length(Dy));
            Dy = Dy/p.delta_xi/2;
            Wellenintervall = 1i*diag(Wellenintervall);
            phi = 1i*phi;
            [Wellenintervall, ind] = sort(Wellenintervall, 'descend', 'ComparisonMethod','real');
            D = diag(Wellenintervall);
            R = phi(:,ind);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            k  = 1;
            dk = 1;
            Q  = 1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif (strcmpi(p.BT,'EXP')==true)
            y = p.xi_Koord;
            Nkk = p.Nk*p.Nk_K;
            Q = zeros(p.N_xi*p.N_K_xi, Nkk);
            k_vec = -2*pi/p.L_xi*((1:(p.Nk*(p.Nk_K-1)+1))-1/2*((p.Nk*(p.Nk_K-1)+2)));
            kDG = linspace(k_vec(1), k_vec(end), p.Nk+1);
            for a = 1 : p.Nk
                k(1, p.Nk_K*(a-1)+1:p.Nk_K*a) = linspace(kDG(a), kDG(a+1), p.Nk_K);
            end    
            dk    = abs(k(2)-k(1));
            
            for IN = 1 : p.Nk*p.Nk_K
                Q(:,IN) = 1./sqrt(p.L_xi)*exp(+1i*k(IN)*y');
            end
            k = diag(k);
            D = k;
            R=Q;
            Dy = 1;
        end


    elseif (p.doFEM == true)
        Dr_xi = Dmatrix1D(p.N_K_xi-1, p.r_xi, p.V_xi);
        Dy = zeros(p.N_xi+1,p.N_xi+1);
        Dy(1:p.N_K_xi,1:p.N_K_xi)=Dr_xi;
        for aa = 1:p.N_xi-1
            Dy(aa+1:aa+p.N_K_xi, aa+1:aa+p.N_K_xi)=[Dy(aa+1,aa+1)+Dr_xi(1,1), Dr_xi(1,p.N_K_xi); 
                                                        Dr_xi(p.N_K_xi,1), Dr_xi(p.N_K_xi,p.N_K_xi)];
        end
        Dy(end,end) = 0;
        Dy(1,1) = 0;
        Dy = p.M_xi_glob\Dy;
        [phi, D] = eig(Dy);
        Wellenintervall = 1i*diag(D);
        phi = 1i*phi;
        Dy = Dy/p.delta_xi/2;
        [D, ind] = sort(Wellenintervall, 'descend', 'ComparisonMethod','real');
        D = diag(D);
        R = phi(:,ind);
    end

end
