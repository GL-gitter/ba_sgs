function [x, w] = JacobiGL(alpha,beta,N)

% function [x, w] = JacobiGL(alpha,beta,N)
% Purpose: Compute the N'th order Gauss Lobatto quadrature
%          points, x, associated with the Jacobi polynomial,
%          of type (alpha,beta) > -1 ( <> -0.5). Additionaly compute
%          corresponding quadrature points w.

x = zeros(N+1,1);
w = zeros(N+1,1);
if (N==0)
    x(1) = 0.0;
    w(1) = 2.0;
    return;
end
if (N==1)
    x(1)=-1.0;
    x(2)=1.0;
    w(1)=1;
    w(2)=1;
    return;
end

[xint,~] = JacobiGQ(alpha+1,beta+1,N-2);
x = [-1, xint', 1]';
n=N+1;
for i=1:ceil(length(x(:))/2)
    w(i) = 2/(n*(n-1)) / (legendreP(n-1, x(i)))^2;
    w(end-(i-1)) = w(i);
end
% assert(abs(sum(w)-2)<1e-15);
return;
