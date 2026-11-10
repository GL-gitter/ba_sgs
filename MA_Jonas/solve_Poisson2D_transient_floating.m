function Vbi_OUT = solve_Poisson2D_transient(mat, Vbi_IN, Ef0, IT, n)

% initializing the data

q = 1.602E-19;
eps0 = 8.854E-12;

Nx = mat.Nx;
Ny = mat.Ny;

dx = mat.dx;
dy = mat.dy;

Epsr=[mat.Eps(1,:); mat.Eps; mat.Eps(Nx,:)];
Epsr=[Epsr(:,1), Epsr, Epsr(:,Ny)];

N = mat.Nd*1E+6;

% assign the node numbering

node_num = zeros(Nx+2, Ny+2);
node_type = zeros(Nx+2, Ny+2);

free_node = zeros(Nx*Ny,1);
bnd_node = zeros(Nx*Ny,1);

node_cnt = 0;

for IX = 2 : Nx+1
    
    for IY = 2 : Ny+1
        
        node_cnt = node_cnt+1;
        
        node_num(IX, IY) = node_cnt;
        node_type(IX, IY) = 0;
        
        if mat.boundary(IX-1, IY-1) 
            
            node_type(IX, IY) = 1;
            
        end
        
    end
    
end

ii = 0;
jj = 0;
v = 0;

n_of_el = 0;

for IX = 2 : Nx+1
    
    for IY = 2 : Ny+1
        
        vC = 0;
        node_act = node_num(IX, IY);
        
        if node_type(IX, IY) == 0
            
            if node_num(IX-1, IY) ~= 0
                
                n_of_el = n_of_el+1;
                
                ii(n_of_el) = node_act;
                jj(n_of_el) = node_num(IX-1, IY);
                
                v(n_of_el) = (eps0/1E-18)*(Epsr(IX-1,IY)+Epsr(IX,IY))/2/dx^2;
                vC = vC-v(n_of_el);
                
            end
            
            % right node
            if node_num(IX+1, IY)~=0
                
                n_of_el = n_of_el+1;
                
                ii(n_of_el) = node_act;
                jj(n_of_el) = node_num(IX+1,IY);
                
                v(n_of_el) = (eps0/1E-18)*(Epsr(IX+1,IY)+Epsr(IX,IY))/2/dx^2;
                vC = vC-v(n_of_el);
                
            end
            
            % upper node
            if node_num(IX, IY-1)~=0
                
                n_of_el = n_of_el+1;
                
                ii(n_of_el) = node_act;
                jj(n_of_el) = node_num(IX,IY-1);
                
                v(n_of_el) = (eps0/1E-18)*(Epsr(IX,IY-1)+Epsr(IX,IY))/2/dy^2;
                vC = vC-v(n_of_el);
                
            end
            
            % lower node
            if node_num(IX, IY+1)~=0
                
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
            
        elseif node_type(IX,IY) == 1
            
            n_of_el = n_of_el+1;
            
            ii(n_of_el) = node_act;
            jj(n_of_el) = node_num(IX, IY);
            
            v(n_of_el) = 1;
            
        end
        
    end
    
end

M = sparse(ii, jj, v, node_cnt, node_cnt);
%figure, spy(M)

node_cnt = 0;

for IX = 1 : Nx
    
    for IY = 1 : Ny
        
        node_cnt = node_cnt+1;
        
        n_p(node_cnt, 1) = n(IX, IY);
        N_p(node_cnt, 1) = N(IX, IY);
        Vbi(node_cnt, 1) = 0;
        
        if mat.boundary(IX, IY)
            
            Vbi(node_cnt, 1) = Ef0 + mat.phi_m_g-mat.Xi-mat.Vg(IT);
            
        end
        
%         if IX == 1
%             
%             Vbi(node_cnt, 1) = Vbi_IN(IX, IY);
%             
%         end
%         
%         if IX == Nx
%             
%             Vbi(node_cnt, 1) = Vbi_IN(IX, IY)-mat.Vd(IT);
%             
%         end
        
    end
    
end

RHS = -q*(N_p-n_p)+M*Vbi;

node_num = node_num(2:Nx+1,2:Ny+1);
node_type = node_type(2:Nx+1, 2:Ny+1);

free_node = node_num(node_type==0);
bnd_node = node_num(node_type==1);

Vbi(free_node) = -M(free_node,free_node)\RHS(free_node);

node_cnt = 0;

for IX = 1 : Nx
    
    for IY = 1 : Ny
        
        node_cnt = node_cnt + 1;
        
        Vbi_OUT(IX, IY) = Vbi(node_cnt);
            
    end
    
end
            
            
            
            
            
        
        

