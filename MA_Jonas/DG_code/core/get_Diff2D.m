function [A,rhs] = get_Diff2D(p)

%% Diffusion matrix
% S_chi = [1/(2*p.dxi2D)-1/(2*p.dchi2D), 1/(2*p.dchi2D)-1/(2*p.dxi2D), 1/(2*p.dxi2D)-1/(2*p.dchi2D), 1/(2*p.dchi2D)-1/(2*p.dxi2D);
%          1/(2*p.dchi2D)-1/(2*p.dxi2D), 1/(2*p.dxi2D)+1/(2*p.dchi2D), 1/(2*p.dchi2D)-1/(2*p.dxi2D), 1/(2*p.dxi2D)-1/(2*p.dchi2D);
%          1/(2*p.dchi2D)+1/(2*p.dxi2D), 1/(2*p.dxi2D)-1/(2*p.dchi2D), 1/(2*p.dchi2D)+1/(2*p.dxi2D), 1/(2*p.dxi2D)-1/(2*p.dchi2D);
%          1/(2*p.dxi2D)-1/(2*p.dchi2D), 1/(2*p.dxi2D)+1/(2*p.dchi2D), 1/(2*p.dxi2D)-1/(2*p.dchi2D), 1/(2*p.dchi2D)-1/(2*p.dxi2D)];

S_chi = [-1/6,-1/6,-1/12,-1/12;
          1/6,1/6,1/12,-1/12;
         -1/12,1/12,1/6,1/6;
         -1/12,-1/12,-1/6,-1/6];
S_chiCell = repmat({S_chi},1,p.N_chi2D);
S_chi = blkdiag(S_chiCell{:});

[~,Fchi2D] = flux2D(p);

% A = S_chi/p.dchi2D+p.D_xi-Fchi2D;

A = kron((diag(p.D) > 0).*(0.5*(p.D+abs(p.D))),S_chi-Fchi2D)...
    +kron((diag(p.D) < 0).*(0.5*(abs(p.D)-p.D)),S_chi-Fchi2D);
A = p.Q_diff/p.dchi2D*sparse(A);

%% Boundaries
detJx = ones(p.N_chi,1)*p.dchi2D/2;
Q1 = zeros(p.N_chi2D*p.N2Dp*p.N2Dp,1);
Q1(1) = 1/detJx(1);
Q1(end) = 1/detJx(end);
Q2 = zeros(p.N_chi2D*p.N2Dp*p.N2Dp,1);
Q2(1) = 1/detJx(1);
Q2(end)= 1/detJx(end);
rhs = kron((diag(p.D) > 0).*(abs(p.D)*p.fermi2D),Q1)...
    + kron((diag(p.D) < 0).*(abs(p.D)*p.fermi2D),Q2);
rhs = p.Q_diff*rhs;

% % get triangle coordinates
% [x,y] = Nodes2D(p.N2Dp-1);
% [r,s] = xytors(x,y);

% % get 2D Vandermonde matrix
% V = Vandermonde2D(p.N2Dp-1,r,s);
% invV = inv(V);

% [Dr,Ds] = Dmatrices2D(p.N2Dp-1,r,s,V);

% Jk = r(1)*s(2)-r(2)*s(1);
% M = Jk*invV.'*invV;

% Sr = M*Dr;
% Ss = M*Ds;

end
