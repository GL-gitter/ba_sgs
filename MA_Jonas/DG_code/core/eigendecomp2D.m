function [D_xi,phi,EW] = eigendecomp2D(p)

    [F_xi2D,~] = flux2D(p);

    S_xiN = [-1/6, -1/12, -1/12, -1/6; 
             -1/12, -1/6, -1/6, -1/12;
              1/12, 1/6, 1/6, 1/12;
             -1/6, -1/12, -1/12, -1/6];
    
    S_xiCell = repmat({S_xiN},1,p.N_xi2D);
    S_xi = blkdiag(S_xiCell{:});

    D_xi = (S_xi-F_xi2D)/p.dxi2D;

    [phi,EW] = eig(D_xi);
    Wellenintervall_a = diag(EW);
%     phi = 1i*phi;
    [Wellenintervall_a,ind] = sort(Wellenintervall_a,'descend','ComparisonMethod','real');
    Wellenintervall = zeros(p.N_xi2D,1);
    for i = 1:p.N_xi2D
        Wellenintervall(i) = sum(Wellenintervall_a((i-1)*4+1:i*4))/4;
    end
    EW = diag(Wellenintervall);
    phi = phi(:,ind);
    
end
