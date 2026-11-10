function [V_xyl, V_xyr] = fevalCoordinateTransformALLMB(x, y, nm, func1)
% function fevalCoordinateTransform returns the rotated function for
% V_xy(r,s) = func1(r+s/2) + func2(r-s/2)
% ATTENTION TO THE PLUS SIGN INCLUDED 
%

Ny = length(y);
Nx = length(x);

V_xyl = zeros(nm,nm, 4, Ny, Nx);
V_xyr = zeros(nm,nm, 4, Ny, Nx);

val1 = zeros(nm,nm, 4);
val2 = zeros(nm,nm, 4);

for IY = 1 : Ny
    
    for IX = 1 : Nx
        
        temp_coord1 = x(IX) + 1/2*y(IY);
        temp_coord2 = x(IX) - 1/2*y(IY);
        
        if temp_coord1 >= x(1) && temp_coord1 <= x(Nx)
            for IMO = 1:nm
                for IMI = 1:nm
                    for IK = 1:4
                        Vleft_bound  = func1(IMO,IMI,IK,temp_coord1>= x);
                        Vright_bound = func1(IMO,IMI,IK,temp_coord1<= x);
                    
                        val1(IMO,IMI,IK) = (Vleft_bound(end) + Vright_bound(1))/2;
                    end
                end
            end

        elseif temp_coord1 <= x(1)
            for IMO = 1:nm
                for IMI = 1:nm
                    for IK = 1:4
                        val1(IMO,IMI,IK) = func1(IMO,IMI,IK,1);
                    end
                end
            end
            
        elseif temp_coord1 >= x(Nx)
            for IMO = 1:nm
                for IMI = 1:nm
                    for IK = 1:4
                        val1(IMO,IMI,IK) = func1(IMO,IMI,IK,Nx);
                    end
                end
            end
        end

        if temp_coord2 >= x(1) && temp_coord2 <= x(Nx)
            for IMO = 1:nm
                for IMI = 1:nm
                    for IK = 1:4
                        Vleft_bound = func1(IMO,IMI,IK,temp_coord2>= x);
                        Vright_bound = func1(IMO,IMI,IK,temp_coord2<= x);
                
                        val2(IMO,IMI,IK) = (Vleft_bound(end) + Vright_bound(1))/2;
                    end
                end
            end
            
        elseif temp_coord2 <= x(1)
            
            for IMO = 1:nm
                for IMI = 1:nm
                    for IK = 1:4
                        val2(IMO,IMI,IK) = func1(IMO,IMI,IK,1);
                    end
                end
            end
            
        elseif temp_coord2 >= x(Nx)
            
            for IMO = 1:nm
                for IMI = 1:nm
                    for IK = 1:4
                        val2(IMO,IMI,IK) = func1(IMO,IMI,IK,Nx);
                    end
                end
            end
            
        end
        
        V_xyl(:,:,:,IY, IX) = val1;
        V_xyr(:,:,:,IY, IX) = val2;

    end
    
end





