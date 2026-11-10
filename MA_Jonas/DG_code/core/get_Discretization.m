function X = get_Discretization(p)

    detJy = ones(p.N_xi,1);
    invMxi = inv(p.M_xi);
    Dr_xi = p.M_xi*Dmatrix1D(p.N_K_xi-1,p.r_xi,p.V_xi);
    Dr_ind = find(abs(Dr_xi)<1e-15);
    Dr_xi(Dr_ind) = 0;

    P1 = zeros(p.N_K_xi);
    P2 = P1;
    P3 = P1;
    P4 = P1;
    P1(1,1) = 0.5;
    P2(1, p.N_K_xi) = 0.5;
    P3(p.N_K_xi,p.N_K_xi) = 0.5;
    P4(p.N_K_xi,1) = 0.5;

    K1 = invMxi*(Dr_xi+0.5*P1-0.5*P3);
    K2 = invMxi*(0.5*P1-0.5*P3);
    K3 = -0.5*invMxi*P2;
    K4 = K3;
    K5 = 0.5*invMxi*P4;
    K6 = K5;

    H1 = kron(spdiags(1./detJy,1,p.N_xi,p.N_xi),K5)+kron(spdiags(1./detJy(2:end),-1,p.N_xi,p.N_xi),K3)...
        +kron(spdiags(1./detJy,0,p.N_xi,p.N_xi),K1);
    H2 = kron(spdiags(1./detJy,1,p.N_xi,p.N_xi),K6)+kron(spdiags(1./detJy(2:end),-1,p.N_xi,p.N_xi),K4)...
        +kron(spdiags(1./detJy,0,p.N_xi,p.N_xi),K2);
    X = H1+H2;
    X_ind = find(abs(X)<1e-15);
    X(X_ind) = 0;
end