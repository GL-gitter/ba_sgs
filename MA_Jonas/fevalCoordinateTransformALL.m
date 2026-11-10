function [V_xyl, V_xyr] = fevalCoordinateTransformALL(x, y, nm, func1)
% function fevalCoordinateTransform returns the rotated function for
% V_xy(r,s) = func1(r+s/2) + func2(r-s/2)
% ATTENTION TO THE PLUS SIGN INCLUDED 
%

Ny = length(y);
Nx = length(x);

V_xyl = zeros(nm, 3, Ny, Nx);
V_xyr = zeros(nm, 3, Ny, Nx);

val1 = zeros(nm, 3);
val2 = zeros(nm, 3);

for IY = 1 : Ny
    
    for IX = 1 : Nx
        
        temp_coord1 = x(IX) + 1/2*y(IY);
        temp_coord2 = x(IX) - 1/2*y(IY);
        
        if temp_coord1 >= x(1) && temp_coord1 <= x(Nx)
            for IM = 1:nm
                for IK = 1:3
                    Vleft_bound  = func1(IM,IK,temp_coord1>= x);
                    Vright_bound = func1(IM,IK,temp_coord1<= x);
                
                    val1(IM,IK) = (Vleft_bound(end) + Vright_bound(1))/2;
                end
            end

        elseif temp_coord1 <= x(1)
            for IM = 1:nm
                for IK = 1:3
                    val1(IM,IK) = func1(IM,IK,1);
                end
            end
            
        elseif temp_coord1 >= x(Nx)
            for IM = 1:nm
                for IK = 1:3
                    val1(IM,IK) = func1(IM,IK,Nx);
                end
            end
        end

        if temp_coord2 >= x(1) && temp_coord2 <= x(Nx)
            for IM = 1:nm
                for IK = 1:3
                    Vleft_bound = func1(IM,IK,temp_coord2>= x);
                    Vright_bound = func1(IM,IK,temp_coord2<= x);
            
                    val2(IM,IK) = (Vleft_bound(end) + Vright_bound(1))/2;
                end
            end
            
        elseif temp_coord2 <= x(1)
            
            for IM = 1:nm
                for IK = 1:3
                    val2(IM,IK) = func1(IM,IK,1);
                end
            end
            
        elseif temp_coord2 >= x(Nx)
            
            for IM = 1:nm
                for IK = 1:3
                    val2(IM,IK) = func1(IM,IK,Nx);
                end
            end
            
        end
        
        V_xyl(:,:,IY, IX) = val1;
        V_xyr(:,:,IY, IX) = val2;

    end
    
end





