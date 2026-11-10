function [Fxi2D,Fchi2D] = flux2D(p)

    % local Flux
    F2Dl = [1/2,0,-1/2,0,0,0,0,0;
            0,-1/2,0,1/2,0,0,0,0;
            0,0,0,0,1/2,0,-1/2,0;
            0,0,0,0,0,-1/2,0,1/2];
    
    Fchi2D = zeros(p.N_chi2D*p.N2Dp*p.N2Dp);
    Fxi2D = zeros(p.N_xi2D*p.N2Dp*p.N2Dp);

    % global Flux
    for i = 1:p.N_chi2D-1
        Fchi2D((i*4)+1:(i*4)+4,(i*4)-1:(i*4)+6) = F2Dl(:,:);
    end 
    for i = 1:p.N_xi2D-1
        Fxi2D((i*4)+1:(i*4)+4,(i*4)-1:(i*4)+6) = F2Dl(:,:);
    end
    
    % fine-tuning Fluxes
    Fchi2D(1:2*p.N2Dp,1:3*p.N2Dp) = F2Dl(:,3:end);
    Fxi2D(1:2*p.N2Dp,1:3*p.N2Dp) = F2Dl(:,3:end);
    for i = 1:2
        Fxi2D(:,end) = [];
        Fchi2D(:,end) = [];
    end

    % inter-elemental coupling
    for i = 1:length(p.Kopplung_x(:,1))
        for j = 1:length(p.Kopplung_x(1,:))/p.N2Dp
            if mod(i,2) ~= 0
                Fchi2D(p.Kopplung_x(i,(j-1)*2+1),p.Kopplung_x(i,(j-1)*2+2)) = F2Dl(3,7);
                Fchi2D(p.Kopplung_x(i,(j-1)*2+2),p.Kopplung_x(i,(j-1)*2+1)) = F2Dl(1,1);
            else
                Fchi2D(p.Kopplung_x(i,(j-1)*2+1),p.Kopplung_x(i,(j-1)*2+2)) = F2Dl(2,2);
                Fchi2D(p.Kopplung_x(i,(j-1)*2+2),p.Kopplung_x(i,(j-1)*2+1)) = F2Dl(end,end);
            end
        end
    end
    for i = 1:length(p.Kopplung_y(:,1))/p.N2Dp
        for j = 1:length(p.Kopplung_y(1,:))
            if mod(j,2) ~= 0
                Fxi2D(p.Kopplung_y((i-1)*2+1,j),p.Kopplung_y((i-1)*2+2,j)) = F2Dl(3,7);
                Fxi2D(p.Kopplung_y((i-1)*2+2,j),p.Kopplung_y((i-1)*2+1,j)) = F2Dl(1,1);
            else
                Fxi2D(p.Kopplung_y((i-1)*2+1,j),p.Kopplung_y((i-1)*2+2,j)) = F2Dl(2,2);
                Fxi2D(p.Kopplung_y((i-1)*2+2,j),p.Kopplung_y((i-1)*2+1,j)) = F2Dl(end,end);
            end
        end
    end
end