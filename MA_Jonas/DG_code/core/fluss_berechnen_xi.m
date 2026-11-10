function [Fluss_Global_pos, Fluss_Global_neg] = fluss_berechnen_xi(N_Knoten, N_chi, alpha)

Fluss_Global_pos = zeros(N_chi*N_Knoten, N_chi*N_Knoten);
Fluss_Global_neg = zeros(N_chi*N_Knoten, N_chi*N_Knoten);

f = 1/2 + (1-alpha)/2;

% Fluss_Lokal_pos = [0, 0, 1,-1, 0, 0, 0, 0, 0;                                            
%                    0, 0, 0, 0, 0, 0, 0, 0, 0; 
%                    0, 0, 0, 0, 0, 0, 0, 0, 0];
               
Fluss_Lokal_pos = zeros(N_Knoten, 3*N_Knoten);
Fluss_Lokal_pos(1,N_Knoten) = -f;
Fluss_Lokal_pos(1,N_Knoten+1) = f;
Fluss_Lokal_pos(N_Knoten, N_Knoten) = 0;
Fluss_Lokal_pos(N_Knoten, N_Knoten+1) = 0;

% Z = zeros(N_Knoten,1);
% Fluss_Lokal_pos = [Z Z invM(:,1) -invM(:,1) Z Z Z Z Z];

for a = 1 : N_chi - 2
    Fluss_Global_pos(N_Knoten*a+1 : N_Knoten*a+N_Knoten , N_Knoten*(a-1)+2 : N_Knoten*(a-1)+1+3*N_Knoten) = Fluss_Lokal_pos;                           % nummerischer Fluss für eine Welle
end

Fluss_Global_pos= Fluss_Global_pos(:,2:N_chi*N_Knoten);
Fluss_Global_pos(1:N_Knoten, 1:N_Knoten)= Fluss_Lokal_pos(1:N_Knoten, N_Knoten+1:N_Knoten*N_Knoten);                                              
Fluss_Global_pos(N_Knoten*N_chi-1 : N_Knoten*N_chi, N_Knoten*N_chi-3 : N_Knoten*N_chi) =...
    Fluss_Lokal_pos(1:N_Knoten, 1:2*N_Knoten);

% Fluss_Lokal_neg = [0, 0, 0, 0, 0, 0, 0, 0, 0;                                            
%                    0, 0, 0, 0, 0, 0, 0, 0, 0; 
%                    0, 0, 0, 0, 0,-1, 1, 0, 0];
               
Fluss_Lokal_neg = zeros(N_Knoten, 3*N_Knoten);
Fluss_Lokal_neg(1,N_Knoten) = 0;
Fluss_Lokal_neg(1,N_Knoten+1) = 0;
Fluss_Lokal_neg(N_Knoten, 2*N_Knoten) = -f;
Fluss_Lokal_neg(N_Knoten, 2*N_Knoten+1) = f;

% Z = zeros(N_Knoten,1);
% Fluss_Lokal_neg = [Z Z Z Z Z -invM(:,3) invM(:,3) Z Z];
           
for a = 1 : N_chi - 2
    Fluss_Global_neg(N_Knoten*a+1 : N_Knoten*a+N_Knoten, N_Knoten*(a-1)+2 : N_Knoten*(a-1)+1+3*N_Knoten) = Fluss_Lokal_neg;                           % nummerischer Fluss für eine Welle
end

Fluss_Global_neg=Fluss_Global_neg(:,2:N_chi*N_Knoten);
Fluss_Global_neg(1:N_Knoten, 1:2*N_Knoten)= Fluss_Lokal_neg(1:N_Knoten, N_Knoten+1:3*N_Knoten);                                              
Fluss_Global_neg(N_Knoten*N_chi-1 : N_Knoten*N_chi, N_Knoten*N_chi-3 : N_Knoten*N_chi) =...
    Fluss_Lokal_neg(1:N_Knoten, 1:2*N_Knoten);
end
