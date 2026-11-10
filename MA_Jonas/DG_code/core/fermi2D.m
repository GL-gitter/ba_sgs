function g = fermi2D(N_xi, phi)
a1 = 4.5097e23;
a2 = 4.0577;
a3 = 1.7861;
L_xi = 50;
delta_xi = L_xi/N_xi;

y = linspace(-(L_xi-delta_xi/2)/2, (L_xi-delta_xi/2)/2, N_xi)';

ffd = @(k) a1*cos(y*k).*log(1+exp(-a2.*k^2+a3));
upperk = 2*sqrt(a3/a2);

g = integral(ffd, -upperk, upperk, 'ArrayValued', true);


g = phi\g;
end

