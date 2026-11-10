function [Fluss_Global_pos, Fluss_Global_neg, Flussxi_Global_pos, Flussxi_Global_neg] = fluss_berechnen(N_Knoten, N_chi, alpha, invM)

Fluss_Global_pos = zeros(N_chi*N_Knoten, N_chi*N_Knoten);
Fluss_Global_neg = zeros(N_chi*N_Knoten, N_chi*N_Knoten);
Flussxi_Global_pos = zeros(N_chi*N_Knoten, N_chi*N_Knoten);
Flussxi_Global_neg = zeros(N_chi*N_Knoten, N_chi*N_Knoten);

f = 1/2 + (1-alpha)/2;

% Fluss_Lokal_pos = [0, 0, 1, -1, 0, 0, 0, 0, 0;                                            
%                    0, 0, 0, 0, 0, 0, 0, 0, 0; 
%                    0, 0, 0, 0, 0, 0, 0, 0, 0];
               
Flussxi_Lokal_pos = zeros(N_Knoten, N_Knoten^2);
Flussxi_Lokal_pos(1,N_Knoten) = f;
Flussxi_Lokal_pos(1,N_Knoten+1) = -f;
Flussxi_Lokal_pos(N_Knoten, 2*N_Knoten) = 0;
Flussxi_Lokal_pos(N_Knoten, 2*N_Knoten+1) = 0;

Z = zeros(N_Knoten,1);
Fluss_Lokal_pos = [Z Z f*invM(:,1) -f*invM(:,1) Z Z Z Z Z];

for a = 1 : N_chi - 2
    Fluss_Global_pos(N_Knoten*a+1 : N_Knoten*a+3 , N_Knoten*(a-1)+1 : N_Knoten*(a-1)+9) = Fluss_Lokal_pos;                           % nummerischer Fluss für eine Welle
    Flussxi_Global_pos(N_Knoten*a+1 : N_Knoten*a+3 , N_Knoten*(a-1)+1 : N_Knoten*(a-1)+9) = Flussxi_Lokal_pos;
end

Fluss_Global_pos(1 : 3, 1:6)= Fluss_Lokal_pos(1:3, 4:9);                                              
Fluss_Global_pos(N_Knoten*N_chi-2 : N_Knoten*N_chi, N_Knoten*N_chi-5 : N_Knoten*N_chi) =...
    Fluss_Lokal_pos(1:3, 1:6);
Flussxi_Global_pos(1 : 3, 1:6)= Flussxi_Lokal_pos(1:3, 4:9);                                              
Flussxi_Global_pos(N_Knoten*N_chi-2 : N_Knoten*N_chi, N_Knoten*N_chi-5 : N_Knoten*N_chi) =...
    Flussxi_Lokal_pos(1:3, 1:6);

% Fluss_Lokal_neg = [0, 0, 0, 0, 0, 0, 0, 0, 0;                                            
%                    0, 0, 0, 0, 0, 0, 0, 0, 0; 
%                    0, 0, 0, 0, 0, 1, -1, 0, 0];
               
Flussxi_Lokal_neg = zeros(N_Knoten, N_Knoten^2);
Flussxi_Lokal_neg(1,N_Knoten) = 0;
Flussxi_Lokal_neg(1,N_Knoten+1) = 0;
Flussxi_Lokal_neg(N_Knoten, 2*N_Knoten) = -f;
Flussxi_Lokal_neg(N_Knoten, 2*N_Knoten+1) = f;

Z = zeros(N_Knoten,1);
Fluss_Lokal_neg = [Z Z Z Z Z -f*invM(:,3) f*invM(:,3) Z Z];
           
for a = 1 : N_chi - 2
    Fluss_Global_neg(N_Knoten*a+1 : N_Knoten*a+3 , N_Knoten*(a-1)+1 : N_Knoten*(a-1)+9) = Fluss_Lokal_neg;                           % nummerischer Fluss für eine Welle
    Flussxi_Global_neg(N_Knoten*a+1 : N_Knoten*a+3 , N_Knoten*(a-1)+1 : N_Knoten*(a-1)+9) = Flussxi_Lokal_neg;
end

Fluss_Global_neg(1 : 3, 1:6)= Fluss_Lokal_neg(1:3, 4:9);                                              
Fluss_Global_neg(N_Knoten*N_chi-2 : N_Knoten*N_chi, N_Knoten*N_chi-5 : N_Knoten*N_chi) =...
    Fluss_Lokal_neg(1:3, 1:6);
Flussxi_Global_neg(1 : 3, 1:6)= Flussxi_Lokal_neg(1:3, 4:9);                                              
Flussxi_Global_neg(N_Knoten*N_chi-2 : N_Knoten*N_chi, N_Knoten*N_chi-5 : N_Knoten*N_chi) =...
    Flussxi_Lokal_neg(1:3, 1:6);

end
