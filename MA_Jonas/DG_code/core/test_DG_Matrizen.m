% Polynomgrad
N = 1;

% Referenzknoten
r = JacobiGL(0,0,N);

% Vandermonde-Matrix
V = Vandermonde1D(N,r);

% Ableitungsmatrix
D = Dmatrix1D(N,r,V);

% Massenmatrix
M = inv(V*V');

% Steifigkeitsmatrix
S = D'*M;

disp('Massenmatrix M:')
disp(M)

disp('Steifigkeitsmatrix S:')
disp(S)