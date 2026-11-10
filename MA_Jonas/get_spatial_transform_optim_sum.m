function B = get_spatial_transform_optim_sum(x, y, V)

B = zeros(length(y), length(x));

for IX = 1 : length(x)
    
    for IY = 1 : length(y)
        
        if x(IX)+1/2*y(IY)>=min(x) && x(IX)+1/2*y(IY)<=max(x)
            
            if ismember(x(IX)+1/2*y(IY), x)
                
                val1 = V((x(IX)+1/2*y(IY))==x);
                
            else
                
                upper_bound = min(find((x(IX)+1/2*y(IY))<x));
                lower_bound = max(find((x(IX)+1/2*y(IY))>x));
                val1 = (V(upper_bound)+V(lower_bound))/2;
                
            end
            
        elseif (x(IX)+1/2*y(IY)) > max(x)
            
            val1 = V(end);
            
        else
            
            val1 = V(1);
            
        end
        if x(IX)-1/2*y(IY)>=min(x) && x(IX)-1/2*y(IY)<=max(x)
            
            if ismember(x(IX)-1/2*y(IY), x)
                
                val2 = V((x(IX)-1/2*y(IY))==x);
                
            else
                
                upper_bound = min(find((x(IX)-1/2*y(IY))<x));
                lower_bound = max(find((x(IX)-1/2*y(IY))>x));
                val2 = (V(upper_bound)+V(lower_bound))/2;
                
            end
            
        elseif (x(IX)-1/2*y(IY)) > max(x)
            
            val2 = V(end);
            
        else
            
            val2 = V(1);
            
        end
        B(IY, IX) = val1+val2;
    end
    
end
