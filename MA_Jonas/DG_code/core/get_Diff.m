function [A,rhs] = get_Diff(p)

%% define flux blending parameter
%%%%%%%%%%%%%%%%%%%%%%%
%alpha=1 is equal to full upwind
%alpha=0 is equal to full central flux
%%%%%%%%%%%%%%%%%%%%%%%
% alpha = round(alpha_profile_interfaces(p.N_chi, 0.005, ceil(0.1*p.N_chi)),4);
alpha = ones((p.N_chi-1),1);
alphaExtL = [1.0;alpha];
alphaExtR = [alpha;1.0];


alphaL_matrix = spdiags(alphaExtL, 0, p.N_chi, p.N_chi);
alphaR_matrix = spdiags(alphaExtR, 0, p.N_chi, p.N_chi);
%%


fnl = p.verteilung_l;
fnr = p.verteilung_r;

Dr = p.M*Dmatrix1D(p.N_K_chi-1, p.r, p.V);
detJx = ones(p.N_chi,1)*p.delta_chi/2;
p.M = p.delta_xi/2*p.M;
%p.M = p.delta_k/2*p.M;
p.invM = inv(p.M);
Dr_index = find(abs(Dr)<1e-15);
Dr(Dr_index) = 0;
P1 = zeros(p.N_K_chi);
P1(1,1)= 1;
P2= zeros(p.N_K_chi);
P2(1,p.N_K_chi)= 1;
P3= zeros(p.N_K_chi);
P3(p.N_K_chi,p.N_K_chi)= 1;
P4= zeros(p.N_K_chi);
P4(p.N_K_chi,1)= 1;

K1= p.invM*(Dr+0.5*P1-0.5*P3);
K2= p.invM*(0.5*P1+0.5*P3);
K2L = p.invM*(0.5*P1);
K2R = p.invM*(0.5*P3);
K3= -0.5*p.invM*P2;
K4= K3;                      %%%%%%% Attention K4 = K3  is default
K5= 0.5*p.invM*P4;
K6= -K5;                     %%%%%%% Attention K6 = -K5 is default


if p.permute_shell == false

    H1= kron(spdiags(1./detJx,1,p.N_chi,p.N_chi),K5)+kron(spdiags(1./detJx(2:end),-1,p.N_chi,p.N_chi),K3)...
        +kron(spdiags(1./detJx,0,p.N_chi,p.N_chi),K1);
    H2= kron(spdiags(1./detJx,1,p.N_chi,p.N_chi),K6)+kron(spdiags(1./detJx(2:end),-1,p.N_chi,p.N_chi),K4)...
        +kron(spdiags(1./detJx,0,p.N_chi,p.N_chi),K2);
    A = p.Q_diff*(kron(p.D,H1)+(1)*kron(abs(p.D),H2));
    
    %% Randbedingungen aufstellen
    
    Q1 = zeros(p.N_K_chi*p.N_chi,1);
    Q1(1:p.N_K_chi) = p.invM(:,1)/detJx(1);
    Q2 = zeros(p.N_K_chi*p.N_chi,1);
    Q2(p.N_K_chi*(p.N_chi-1)+1:p.N_K_chi*p.N_chi)= p.invM(:,p.N_K_chi)/detJx(end);
    rhs = kron((diag(p.D) > 0).*(abs(p.D)*fnl),Q1)...   %abs(p.D)*
        +kron((diag(p.D) < 0).*(abs(p.D)*fnr),Q2);      %abs(p.D)*
    rhs = p.Q_diff*rhs;

elseif p.permute_shell == true
    H1= kron(spdiags(1./detJx,1,p.N_chi,p.N_chi),K5)+kron(spdiags(1./detJx(2:end),-1,p.N_chi,p.N_chi),K3)...
        +kron(spdiags(1./detJx,0,p.N_chi,p.N_chi),K1);
    % H2= kron(spdiags(1./detJx,1,p.N_chi,p.N_chi),K6)+kron(spdiags(1./detJx(2:end),-1,p.N_chi,p.N_chi),K4)...
    %     +kron(spdiags(1./detJx,0,p.N_chi,p.N_chi),K2);
    H2= kron(spdiags(1./detJx,1,p.N_chi,p.N_chi)*alphaL_matrix,K6)+kron(spdiags(1./detJx(2:end),-1,p.N_chi,p.N_chi)*alphaR_matrix,K4)...
         +kron(spdiags(1./detJx,0,p.N_chi,p.N_chi)*alphaL_matrix,K2L)+kron(spdiags((1./detJx),0,p.N_chi,p.N_chi)*alphaR_matrix,K2R);
    A = p.Q_diff*(kron(H1,p.D)+kron(H2,abs(p.D)));
    
    %% Randbedingungen aufstellen
    
    Q1 = zeros(p.N_K_chi*p.N_chi,1);
    Q1(1:p.N_K_chi) = p.invM(:,1)/detJx(1); 
    Q2 = zeros(p.N_K_chi*p.N_chi,1);
    Q2(p.N_K_chi*(p.N_chi-1)+1:p.N_K_chi*p.N_chi)= p.invM(:,p.N_K_chi)/detJx(end);
    rhs = kron(Q1,(diag(p.D) > 0).*(abs(p.D)*fnl))...   %abs(p.D)*
        +kron(Q2,(diag(p.D) < 0).*(abs(p.D)*fnr));      %abs(p.D)*
    rhs = p.Q_diff*rhs;


end
end






%% former version



%XX = diag(max(max(abs(p.D)))*ones(p.Nk,1));     %%%%just try this!

% H1= kron(spdiags(1./detJx,1,p.N_chi,p.N_chi),K5)+kron(spdiags(1./detJx(2:end),-1,p.N_chi,p.N_chi),K3)...
%     +kron(spdiags(1./detJx,0,p.N_chi,p.N_chi),K1);
% H2= kron(spdiags(1./detJx,1,p.N_chi,p.N_chi),K6)+kron(spdiags(1./detJx(2:end),-1,p.N_chi,p.N_chi),K4)...
%     +kron(spdiags(1./detJx,0,p.N_chi,p.N_chi),K2);
% A = p.Q_diff*(kron(p.D,H1)+kron(abs(p.D),H2));
% 
% %% Randbedingungen aufstellen
% 
% Q1 = zeros(p.N_K_chi*p.N_chi,1);
% Q1(1:p.N_K_chi) = p.invM(:,1)/detJx(1);
% Q2 = zeros(p.N_K_chi*p.N_chi,1);
% Q2(p.N_K_chi*(p.N_chi-1)+1:p.N_K_chi*p.N_chi)= p.invM(:,p.N_K_chi)/detJx(end);
% rhs = kron((diag(p.D) > 0).*(abs(p.D)*fnl),Q1)...   %abs(p.D)*
%     +kron((diag(p.D) < 0).*(abs(p.D)*fnr),Q2);      %abs(p.D)*
% rhs = p.Q_diff*rhs;
% 
% end


%diag(max(max(abs(p.D)))*ones(p.Nk,1));