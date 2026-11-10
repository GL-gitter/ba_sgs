function [Vp] = V_transient(x, xg, Ug, U, Lx)

Vp = zeros(size(x)); 
Vp(x < Lx)= interp1(xg, Ug, x(x < Lx));
Vp(x >= Lx)= min(Ug);
Vp(x < 0)= Ug(1);

end
