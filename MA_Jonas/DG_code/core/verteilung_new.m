function g = verteilung(p, Ny)
    % a1 = 4.5097e23;
    % a2 = 4.0577;
    % a3 = 1.7861;
    % L_xi = 50;
    % delta_xi = L_xi/N_xi;
    % 
    % y = linspace(-(L_xi-delta_xi/2)/2, (L_xi-delta_xi/2)/2, N_xi)';main
    % 
    % ffd = @(k) a1*cos(y*k).*log(1+exp(-a2.*k^2+a3));
    % upperk = 2*sqrt(a3/a2);
    % 
    % g = integral(ffd, -upperk, upperk, 'ArrayValued', true);
    % 
    % g = phi\g;

    EfL = 1.5090;
    EfR = 1.5090;

    k_fermi = -2*pi/p.L_xi*((1:Ny)-1/2*(Ny+1));
    dk_fermi = abs(k_fermi(2)-k_fermi(1));

    EL = 1.4240+p.hquer^2.*k_fermi.^2./2/p.m_elektron/p.q;
    ER = 1.4240+p.hquer^2*k_fermi.^2/2/p.m_elektron/p.q;
    fL = 2*p.m_elektron*p.kB*p.T/(2*pi*p.hquer^2)*log(1+exp(p.q*(EfL-EL)/p.kB/p.T));
    fR = 2*p.m_elektron*p.kB*p.T/(2*pi*p.hquer^2)*log(1+exp(p.q*(EfR-ER)/p.kB/p.T));

    fLn = cos(kron(p.xi_Koord',k_fermi))*fL.'*dk_fermi/2/pi;
    fRn = cos(kron(p.xi_Koord',k_fermi))*fR.'*dk_fermi/2/pi;

    fLn = p.phi'*fLn;
    fRn = p.phi'*fRn;
    g = [fLn;fRn];
end

