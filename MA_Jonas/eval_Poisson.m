function [M, node_num, node_type] = eval_Poisson(mat)

q = 1.602E-19;
eps0 = 8.854E-12;

Nx = mat.Nx;
Ny = mat.Ny;

dx = mat.dx;
dy = mat.dy;

Epsr=[mat.Eps(1,:); mat.Eps; mat.Eps(Nx,:)];
Epsr=[Epsr(:,1), Epsr, Epsr(:,Ny)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assign the node numbering

node_num  = zeros(Nx+2, Ny+2);
node_type = zeros(Nx+2, Ny+2); 
node_cnt = 0;

for IX = 2 : Nx+1
    
    for IY = 2 : Ny+1
        
        node_cnt = node_cnt+1;
        node_num(IX, IY) = node_cnt;
        node_type(IX, IY) = 1+mat.boundary(IX-1,IY-1);
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assembling of the node coupling matrix

ii = 0;
jj = 0;
v = 0;
n_of_el = 0;


for IX = 2 : Nx+1
    
    for IY = 2 : Ny+1
        
        node_act = node_num(IX, IY);
        vC       = 0;
        
        if node_type(IX,IY)==1
            
            % left node 
            if node_type(IX-1, IY)~=0
                
                n_of_el = n_of_el+1;
                ii(n_of_el) = node_act;
                jj(n_of_el) = node_num(IX-1,IY);
                v(n_of_el) = (eps0/1E-18)*(Epsr(IX-1,IY)+Epsr(IX,IY))/2/dx^2;
                vC = vC-v(n_of_el);
                
            end
            
            % right node 
            if node_type(IX+1, IY)~=0
                
                n_of_el = n_of_el+1;
                ii(n_of_el) = node_act;
                jj(n_of_el) = node_num(IX+1,IY);
                v(n_of_el) = (eps0/1E-18)*(Epsr(IX+1,IY)+Epsr(IX,IY))/2/dx^2;
                vC = vC-v(n_of_el);
                
            end
            
            % upper node 
            if node_type(IX, IY-1)~=0
                
                n_of_el = n_of_el+1;
                ii(n_of_el) = node_act;
                jj(n_of_el) = node_num(IX,IY-1);
                v(n_of_el) = (eps0/1E-18)*(Epsr(IX,IY-1)+Epsr(IX,IY))/2/dy^2;
                vC = vC-v(n_of_el);
                
            end
            
            % lower node 
            if node_type(IX, IY+1)~=0
                
                n_of_el = n_of_el+1; % increase number of elems
                ii(n_of_el) = node_act;
                jj(n_of_el) = node_num(IX,IY+1); % coupling node
                v(n_of_el) = (eps0/1E-18)*(Epsr(IX,IY+1)+Epsr(IX,IY))/2/dy^2;
                vC = vC-v(n_of_el); % direct assignment of neumann bc's
                
            end
            
            % central node
            n_of_el = n_of_el+1;
            ii(n_of_el) = node_act;
            jj(n_of_el) = node_num(IX, IY);
            v(n_of_el) = vC;
            
        elseif node_type(IX,IY)==2
            
            n_of_el = n_of_el+1;
            ii(n_of_el) = node_act;
            jj(n_of_el) = node_num(IX, IY);
            v(n_of_el) = 1;
            
        end
        
    end
    
end
        
M = sparse(ii, jj, v, node_cnt, node_cnt);

