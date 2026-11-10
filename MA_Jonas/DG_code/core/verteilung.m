function [fnl,fnr] = verteilung(mat, p, N_xi, phi, EfL, EfR, V, m)

%% general parameters
Nx = mat.Nx;

h = 6.626E-34;
m0 = 9.109E-31;
hb = h/2/pi;
kB = 1.38E-23;
q = 1.602E-19;
Temp = mat.Temp;
%%

if (strcmpi(p.choose_fermi,'MJ')==true)
%% Nach Matthias Jaeger
 a1 = 4.5097e23;
 a2 = 4.0577;
 a3 = 1.7861;
 L_xi = 50;
 delta_xi = L_xi/N_xi;
 
 y = linspace(-(L_xi-delta_xi/2)/2, (L_xi-delta_xi/2)/2, N_xi)';
 
 ffd = @(k) a1*cos(y*k).*log(1+exp(-a2.*k^2+a3));
 upperk = 2*sqrt(a3/a2);
 
 g = integral(ffd, -upperk, upperk, 'ArrayValued', true);
 
 fnl = phi\g;
 fnr = 0;
elseif (strcmpi(p.choose_fermi,'A')==true)
%% Nach Asena
 beta = 1/(p.kB*p.T);
 tau = p.q/p.hquer;
 gamma = sqrt(tau*p.m_elektron/p.hquer);
 a1= 2*gamma*p.m_elektron/(2*pi*p.hquer)^2/beta;
 a2= p.q*beta/2;
 a3= 0.0467*p.q*beta;
 a3=1.5090                  %%%%%%%%%%%%%
 
 k       = -2*pi/p.L_xi*1e-9*((1:N_xi)-1/2*(N_xi+1));
 dk      = abs(k(2)-k(1));
 fk      = a1*log(1+exp(-a2.*k.^2+a3));
 fLn     = dk*cos(kron(k,p.xi_Rechengebiet_Elem'*1e9))*fk';
 fnl      = phi\fLn;
 fnr     = 0;

elseif (strcmpi(p.choose_fermi,'LS')==true)
%% Nach Lukas Schulz
 Ny = N_xi+1-1;
 if (p.doDG == true)
     y = p.xi_Koord;
 else
     y = p.xi_Rechengebiet_Elem;
 end
 EfL = 1.509; %1.7465;
 Vex = EfL-0.0465; %1.7;

% Berechnen der Fermi-Welle-Vektoren
 k_fermi = -2*pi/p.L_xi*((1:N_xi)-1/2*(N_xi+1));
 dk = abs(k_fermi(2)-k_fermi(1));

% Bestimmen der Differenz der Fermi-Welle-Vektoren
 EL   = Vex+p.hquer^2.*k_fermi.^2./2/(p.m_elektron)/p.q;
 fL = 2*p.m_elektron*p.kB*p.T/(2*pi*p.hquer^2)*log(1+exp(p.q*(EfL-EL)/p.kB/p.T));

% Transformieren der Fermi-Verteilungen vom Phasenraum in den Ortsraum
 fLn = zeros(N_xi,1);
 for IY = 1:Ny-1
     fLn(IY,1) = cos(y(IY)*k_fermi)*fL.'*dk/2/pi;
 end
 fnl = phi\fLn;
 fnr = 0;


%% Fermi for 3D
elseif (strcmpi(p.choose_fermi,'DG_MSA_3D')==true)
 Ny = N_xi+1-1;
 if (p.doDG == true)
     y = p.xi_Koord;
 else
     y = p.xi_Rechengebiet_Elem;
 end
 %EfL = 1.509; %1.7465;
 %Vex = EfL-0.0465; %1.7;

% Berechnen der Fermi-Welle-Vektoren
 k_fermi = -2*pi/p.L_xi*((1:p.Nk)-1/2*(p.Nk+1));
 dk = abs(k_fermi(2)-k_fermi(1));

%%
 EL   = V(1)+hb^2.*k_fermi.^2./2/m(1)/q;
 ER   = V(Nx)+hb^2*k_fermi.^2/2/m(Nx)/q;

 fL = zeros(1, p.Nk);
 fR = zeros(1, p.Nk);
        
 fL = 2*p.deg_factor_temp*(1./(1+exp(q*(EL-EfL)/(kB*Temp))));
 fR = 2*p.deg_factor_temp*(1./(1+exp(q*(ER-EfR)/(kB*Temp))));
%%

% Bestimmen der Differenz der Fermi-Welle-Vektoren
 %EL   = Vex+p.hquer^2.*k_fermi.^2./2/(p.m_elektron)/p.q;
 %fL   = 2*p.m_elektron*p.kB*p.T/(2*pi*p.hquer^2)*log(1+exp(p.q*(EfL-EL)/p.kB/p.T));

% Transformieren der Fermi-Verteilungen vom Phasenraum in den Ortsraum
 fLn = zeros(N_xi,1);
 fRn = zeros(N_xi,1);
 for IY = 1:Ny
     fLn(IY,1) = cos(y(IY)*k_fermi)*fL.'*dk/2/pi;
     fRn(IY,1) = cos(y(IY)*k_fermi)*fR.'*dk/2/pi;
 end
 fnl = phi\fLn;
 fnr = phi\fRn;
%%

%% Fermi for 2D
elseif (strcmpi(p.choose_fermi,'DG_MSA_2D')==true)
 Ny = N_xi+1-1;
 if (p.doDG == true)
     y = p.xi_Koord;
 else
     y = p.xi_Rechengebiet_Elem;
 end

 if ((strcmpi(p.BT,'EXP') && p.doDG)==true)
   k_fermi = diag(p.k).';     %-2*pi/p.L_xi*((1:p.Nk)-1/2*(p.Nk+1));      %% limit k values here???
   dk = abs(k_fermi(2)-k_fermi(1));
 else
   k_fermi = -2*pi/p.L_xi*((1:p.Nk)-1/2*(p.Nk+1));    %% limit k values here???
   dk = abs(k_fermi(2)-k_fermi(1));
 end

 me_z = mat.me_z_ch;
 EL   = V(1)+hb^2.*k_fermi.^2./2/m(1)/q;
 ER   = V(Nx)+hb^2*k_fermi.^2/2/m(Nx)/q;

  
 kz_max = sqrt(2*2*m0*me_z*q*(max(V(:)+0.3)))/hb;
 
 Nkz = 1001;
 kz = linspace(0, kz_max, Nkz);
 dkz = kz(2)-kz(1);
 Ekz = hb^2*kz.^2/2/m0/me_z/q;
 
 fL = zeros(1,p.Nk*p.Nk_K); 
 fR = zeros(1,p.Nk*p.Nk_K); 

 for IK =1:p.Nk*p.Nk_K
    fL(IK) = 2/pi*dkz*sum(1./(1+exp(q*(EL(IK)+Ekz-EfL)/(kB*Temp))));
    fR(IK) = 2/pi*dkz*sum(1./(1+exp(q*(ER(IK)+Ekz-EfR)/(kB*Temp))));
 end
 
 % Transformieren der Fermi-Verteilungen vom Phasenraum in den Ortsraum
 fLn = zeros(N_xi,1);
 fRn = zeros(N_xi,1);

 for IY = 1:Ny
    fLn(IY,1) = cos(y(IY)*k_fermi)*fL.'*dk/2/pi;
    fRn(IY,1) = cos(y(IY)*k_fermi)*fR.'*dk/2/pi;
 end

 if ((strcmpi(p.BT,'EXP') && p.doDG)==true)
     fnl = phi'*fLn;
     fnr = phi'*fRn;
 else
     fnl = phi\fLn;
     fnr = phi\fRn;
 end
 %%  
end
end
