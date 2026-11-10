function Ht = eval_transversal_Hamiltonian(valley_index, x_index, mat)

h = 6.626E-34;
hb = h/2/pi;
m0 = 9.109E-31;
q = 1.602E-19;
Ny = mat.Ny;
dy = mat.dy;

m_eff = squeeze(mat.me_y(valley_index, x_index, :));

% m_eff=m_eff*0+0.023

V     = mat.V(x_index, :);

% m_eff(6) = m_eff(7);
% m_eff(Ny-5) = m_eff(Ny-6);
% 
% V(6) = V(7);
% V(Ny-5) = V(Ny-6);

t0 = hb^2/2/(m0*dy^2*q*1E-18);

Ht = sparse(Ny, Ny);

for IY = 2 : Ny-1
    
    Ht(IY, IY-1) = -2*t0/(m_eff(IY-1)+m_eff(IY));
    Ht(IY, IY+1) = -2*t0/(m_eff(IY+1)+m_eff(IY));
    
    Ht(IY, IY) = -1*(Ht(IY, IY-1)+Ht(IY,IY+1))+V(IY);
    
end

Ht(1,1) = Ht(2,2); Ht(1,2) = Ht(2,3);
Ht(Ny,Ny) = Ht(Ny-1,Ny-1);
Ht(Ny, Ny-1) = Ht(Ny-1, Ny-2);

%x12345=1;


%plot(kf,1.424+hb^2*kf.^2/2/me/q-hb^2*kf.^4/2/(me*4e18)/q+hb^2*kf.^6/2/(me*5e37)/q)



