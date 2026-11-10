function [V_xy, a, b] = Schwerpunkttrafo(n,func1)
% function fevalCoordinateTransform returns the rotated function for
% V_xy(r,s) = func1(r+s/2) + func2(r-s/2)
% ATTENTION TO THE PLUS SIGN INCLUDED 
%
V_xy = zeros(n,n);

for x = 1:n
    for c = 1:n
        r = (c-1)*1/2+1+1/2*(x-1);
        r_strich = n/2+1/2-1/2*(c-1)+1/2*(x-1);
        if (mod(r,1) == 0 && mod(r_strich,1) == 0) 
            temp = func1(r_strich,r);
        else         
            temp = 1/4*(func1(r_strich-0.5,r-0.5)+func1(r_strich-0.5,r+0.5)+func1(r_strich+0.5,r-0.5)+func1(r_strich+0.5,r+0.5));
            
        end
        V_xy(x,c) = temp;
        a(x,c) = r;
        b(x,c) = r_strich;
    end
end

end