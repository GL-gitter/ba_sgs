function Ht = eval_transversal_Hamiltonian_3D(valley_index, x_index, mat)

h = 6.626E-34;
hb = h/2/pi;
m0 = 9.109E-31;
q = 1.602E-19;
Ny = mat.Ny;
dy = mat.dy;
Nz = mat.Nz;
dz = mat.dz;
m_y = squeeze(mat.me_y(valley_index, x_index, :, :));
m_z = squeeze(mat.me_z(valley_index, x_index, :, :));
% m_eff=m_eff*0+0.023

V     = squeeze(mat.V(x_index, :, :));


Dy2 = zeros(Ny*Nz, Ny*Nz);
Dz2 = zeros(Ny*Nz, Ny*Nz);

Hty = zeros(Ny, Ny);
Htz = zeros(Ny, Ny);
Htz_u = zeros(Ny, Ny);
Htz_l = zeros(Ny, Ny);

t0y = hb^2 / (2 * m0 * dy^2 * q * 1E-18);
t0z = hb^2 / (2 * m0 * dz^2 * q * 1E-18);


% DY2
for IZ = 1:Nz
    for IY = 2:Ny-1
        Hty(IY, IY-1) = -2 * t0y / (m_y(IY-1, IZ) + m_y(IY, IZ));
        Hty(IY, IY+1) = -2 * t0y / (m_y(IY+1, IZ) + m_y(IY, IZ));
        Hty(IY, IY) = -1 * (Hty(IY, IY-1) + Hty(IY, IY+1)) + V(IY, IZ);
    end
    Hty(1, 1) = Hty(2, 2);
    Hty(1, 2) = Hty(2, 3);
    Hty(Ny, Ny) = Hty(Ny-1, Ny-1);
    Hty(Ny, Ny-1) = Hty(Ny-1, Ny-2);
    Dy2((IZ-1)*Ny+1:IZ*Ny, (IZ-1)*Ny+1:IZ*Ny) = Hty;
end

% DZ2
for IZ = 2:Nz-1
    for IY = 1:Ny
        Htz_u(IY, IY) = -2 * t0z / (m_z(IY, IZ-1) + m_z(IY, IZ));
        Htz_l(IY, IY) = -2 * t0z / (m_z(IY, IZ+1) + m_z(IY, IZ));
        Htz(IY, IY) = -1 * (Htz_u(IY, IY) + Htz_l(IY, IY));
    end
    Htz(1, 1) = Htz(2, 2);
    Htz(Ny, Ny) = Htz(Ny-1, Ny-1);
    Htz_l(1, 1) = Htz_l(2, 2);
    Htz_l(Ny, Ny) = Htz_l(Ny-1, Ny-1);
    Htz_u(1, 1) = Htz_u(2, 2);
    Htz_u(Ny, Ny) = Htz_u(Ny-1, Ny-1);
    Dz2((IZ-1)*Ny+1:IZ*Ny, (IZ-1)*Ny+1:IZ*Ny) = Htz;
    Dz2((IZ-1)*Ny+1:IZ*Ny, (IZ-2)*Ny+1:(IZ-1)*Ny) = Htz_u;
    Dz2((IZ-1)*Ny+1:IZ*Ny, IZ*Ny+1:(IZ+1)*Ny) = Htz_l;
end

Dz2(1:Ny, 1:Ny) = Dz2(Ny+1:2*Ny, Ny+1:2*Ny);
Dz2(1:Ny, Ny+1:2*Ny) = Dz2(Ny+1:2*Ny, 2*Ny+1:3*Ny);
Dz2((Nz-1)*Ny+1:Nz*Ny, (Nz-1)*Ny+1:Nz*Ny) = Dz2((Nz-2)*Ny+1:(Nz-1)*Ny, (Nz-2)*Ny+1:(Nz-1)*Ny);
Dz2((Nz-1)*Ny+1:Nz*Ny, (Nz-2)*Ny+1:(Nz-1)*Ny) = Dz2((Nz-2)*Ny+1:(Nz-1)*Ny, (Nz-3)*Ny+1:(Nz-2)*Ny);

Ht = Dy2 + Dz2;

if max(max(abs(Ht - transpose(conj(Ht))))) > 1e-14
    fprintf('Ht not hermitian');
end

% 
% % m_eff(6) = m_eff(7);
% % m_eff(Ny-5) = m_eff(Ny-6);
% % 
% % V(6) = V(7);
% % V(Ny-5) = V(Ny-6);
% 
% t0 = hb^2/2/(m0*dy^2*q*1E-18);
% 
% Ht = sparse(Ny, Ny);
% 
% for IY = 2 : Ny-1
% 
%     Ht(IY, IY-1) = -2*t0/(m_eff(IY-1)+m_eff(IY));
%     Ht(IY, IY+1) = -2*t0/(m_eff(IY+1)+m_eff(IY));
% 
%     Ht(IY, IY) = -1*(Ht(IY, IY-1)+Ht(IY,IY+1))+V(IY);
% 
% end
% 
% Ht(1,1) = Ht(2,2); Ht(1,2) = Ht(2,3);
% Ht(Ny,Ny) = Ht(Ny-1,Ny-1);
% Ht(Ny, Ny-1) = Ht(Ny-1, Ny-2);

%x12345=1;


%plot(kf,1.424+hb^2*kf.^2/2/me/q-hb^2*kf.^4/2/(me*4e18)/q+hb^2*kf.^6/2/(me*5e37)/q)



