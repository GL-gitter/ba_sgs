function CAP = fevalComplexAbsorbingPotential(y)
% evaluates the real space complex absorbing potential

Ny = length(y);
CAP = zeros(1, Ny);
y_max = max(y);

b_al = 4;
w_al = y_max*0.25;
n    = 4;

for i_y = 1 : Ny
    
    if abs(y(i_y)) > y_max - w_al
        
        CAP(i_y) = b_al*((abs(y(i_y))-(y_max - w_al))/w_al)^n;
        
    else
        
        CAP(i_y) = 0;
        
    end
    
end
    
    

