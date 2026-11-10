function [x,n] = getN(p,rho,N)

s= linspace(-1,1,N);


x= linspace(0,p.L_chi,(N-1)*p.N_chi+1);
n= zeros(1,(N-1)*p.N_chi+1);


Vs= Vandermonde1D(p.N_K_chi-1,s);
for k=1:p.N_chi
    rhok= 0.5*(rho((p.N_K_chi*(k-1)+1):(p.N_K_chi*k),p.N_xi*p.N_K_xi/2)+rho((p.N_K_chi*(k-1)+1):(p.N_K_chi*k),p.N_xi*p.N_K_xi/2+1));
    rhok= real(Vs*p.invV*rhok(:));
    n(((k-1)*N+1-(k-1)):(k*N-(k-1)))= rhok;
end

end
