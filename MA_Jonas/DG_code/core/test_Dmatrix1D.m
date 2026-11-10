% Zwei Knoten, also Polynomgrad 1
N = 1;

% Knoten des Referenzelements
r = JacobiGL(0,0,N);

% Vandermonde-Matrix
V = Vandermonde1D(N,r);

% Ableitungsmatrix
D = Dmatrix1D(N,r,V);

% Beispiel: u(r) = 2r + 3
u = 2*r + 3;

% Numerische Ableitung an beiden Knoten
du = D*u;

disp('Knotenpositionen:')
disp(r)

disp('Ableitungsmatrix:')
disp(D)

disp('Berechnete Ableitungen:')
disp(du)