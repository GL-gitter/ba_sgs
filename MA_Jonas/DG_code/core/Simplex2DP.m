function [P] = Simplex2DP(r,s,i,j)

h1 = JacobiP(r,0,0,i);
h2 = JacobiP(s,0,0,j);
P = h1.*h2;
return;
