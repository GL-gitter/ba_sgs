% STUFE 4: Energieaufgeloeste Dichtematrix mit DG(X) + FV(xi)
% ERSTER CUT-CELL-PROTOTYP – noch nicht gegen eine MATLAB-Ausfuehrung validiert.
%
% Voraussetzung: Schritt3.m ausfuehren (phi_werte, E_ABS_plus_num,
% E_ABS_minus_num, M_L, M_R, L, vF, sigma_x, hbar, U0, mu).
% Schritt2.m ist fuer dieses Skript nicht erforderlich.
%
% Ziel: fuer jedes vorgegebene E_ABS(i,n) die beiden SPEKTRALGLEICHUNGEN
% auf der Raute 0<=x1,x2<=L diskretisieren und einen Vektor u finden,
% der beide Gleichungen und die kohärenten Andreev-Randfluesse erfuellt.
%
% WICHTIG: Vor dem gesamten Phasensweep erst einen Index testen!

%% 1. Gitter: DG Grad 1 in X, FV in xi; alle Koordinaten dimensionslos
N = 1;
assert(N==1, 'Dieser Prototyp ist fuer DG-Polynomgrad N=1 ausgelegt.');
K = 8;
N_xi = 17;

% Fuer einen Test nur einen Phasenindex auswaehlen:
i_werte = 175;
% Fuer alle Phasen (sehr rechenaufwendig!) stattdessen:
% i_werte = 1:numel(phi_werte);

E_ABS = [E_ABS_plus_num(:), E_ABS_minus_num(:)];
rho_E_alle = cell(2,numel(phi_werte));
u_alle = cell(2,numel(phi_werte));
residuen = nan(2,numel(phi_werte));

%% 2. Cut-Cell-Geometrie EINMAL erzeugen
% t=X/L, z=xi/L; die Raute ist |z|<=2 min(t,1-t)
[Zellen,idmap] = rauten_gitter(K,N_xi);
N_zellen = numel(Zellen);
N_u = 16*2*N_zellen;
fprintf('Physikalisch aktive Cut-Cells: %d; Unbekannte: %d\n',N_zellen,N_u);

%% 3. Bekannte BdG- und Transportmatrizen
I2 = eye(2);
I4 = eye(4);
tau_z = diag([1,-1]);

A0 = kron(tau_z,sigma_x);       % dimensionslos: A/vF
B = kron(tau_z,(U0-mu)*I2);   % eV, ky=0, konstantes Potential
E_skala = hbar*vF/L;           % eV, damit der Operator dimensionslos wird

% Eigenwertprojektionen von A0 (bei ky=0: Eigenwerte +/-1)
Pplus = (I4+A0)/2;
Pminus = (I4-A0)/2;

% Linkswirkender und rechtswirkender Transport, jeweils 16x16
C1 = kron(I4,A0);
C1p = kron(I4,Pplus);
C1m = kron(I4,-Pminus);       % Negativanteil von A0

C2 = kron(A0.',I4);
C2p = kron(Pplus.',I4);
C2m = kron(-Pminus.',I4);

% Projektionen auf aus-/einlaufende Charakteristiken
P1plus = kron(I4,Pplus);
P1minus = kron(I4,Pminus);
P2plus = kron(Pplus.',I4);
P2minus = kron(Pminus.',I4);

%% 4. Schleife ueber PHASEN und die BEIDEN ABS-ZWEIGE
for i = i_werte
    phi = phi_werte(i);

    for n = 1:2
        E_n = E_ABS(i,n);
        ML = M_L(E_n,phi);
        MR = M_R(E_n,phi);

        % Randbedingungen als 8x16-Matrizen
        R1L = kron(I4,[-ML,I2]);
        R1R = kron(I4,[-MR,I2]);
        Q_L = [-ML'; I2];
        Q_R = [-MR'; I2];
        R2L = kron(Q_L.',I4);
        R2R = kron(Q_R.',I4);

        % Abbildung des einlaufenden Anteils auf auslaufende Randspuren.
        % An x1=0 bzw. x2=0 ist die +Charakteristik einlaufend;
        % an x1=L bzw. x2=L die -Charakteristik.
        T1L = rand_abbildung(R1L,P1plus,P1minus);
        T1R = rand_abbildung(R1R,P1minus,P1plus);
        T2L = rand_abbildung(R2L,P2plus,P2minus);
        T2R = rand_abbildung(R2R,P2minus,P2plus);

        Rand.T1L = T1L; Rand.T1R = T1R;
        Rand.T2L = T2L; Rand.T2R = T2R;

        % Potentialterme beider Spektralgleichungen, vec-Konvention
        P_E = (B-E_n*I4)/E_skala;
        V1 = kron(I4,P_E);       % vec((B-E)*rho)
        V2 = kron(P_E.',I4);     % vec(rho*(B-E))

        %% 5. Beide Spektralmatrizen ueber ALLE Zellen montieren
        [K_Spektral_1,K_Spektral_2] = montiere_cutcell( ...
            Zellen,idmap,C1,C1p,C1m,C2,C2p,C2m,V1,V2,Rand);

        % Beide Gleichungen wirken auf genau denselben Vektor u
        K_gesamt = [K_Spektral_1;K_Spektral_2];

        %% 6. Kleinstes singulaeres Tripel -> Kandidat fuer u
        % Eine DISKRETE ABS-Energie muss kein exakter Nullwert sein.
        try
            [~,Dmin,Vmin] = svds(K_gesamt,1,'smallest');
            u = Vmin(:,1);
            sigma_min = abs(Dmin(1,1));
        catch ME
            warning('svds nicht konvergiert (%s); Versuch ueber K''*K.',ME.message);
            [Vmin,eval_min] = eigs(K_gesamt'*K_gesamt,1,'smallestabs');
            u = Vmin(:,1);
            sigma_min = sqrt(max(0,real(eval_min(1,1))));
        end

        residuum = norm(K_gesamt*u) / max(normest(K_gesamt)*norm(u),eps);
        residuen(n,i) = residuum;
        fprintf('phi-Index %d, Zweig %d: E=%.5g eV, sigma_min=%.3e, rel. Residuum=%.3e\n', ...
            i,n,E_n,sigma_min,residuum);

        %% 7. Die physikalische Dichte normieren: Integral tr(rho(x,x)) dx = 1
        u = normiere_rho(u,Zellen,idmap,K,N_xi,L);

        %% 8. Ueber alle Zellknoten 4x4-Matrizen zurueckgewinnen
        rho_E = reshape(u,4,4,2,N_zellen);
        rho_E_alle{n,i} = rho_E;
        u_alle{n,i} = u;

        % Achtung: Die physikalische Loesung muss per Vergleich mit der
        % aus Schritt4.m bekannten Referenz rho=Psi*Psi' validiert werden.
    end
end

%% 9. Ein Beispiel darstellen: alle 16 Komponenten auf der Raute
% (Nur sofern die ausgewaehlte Berechnung gelungen ist.)
i_plot = i_werte(1);
n_plot = 1;
if ~isempty(rho_E_alle{n_plot,i_plot})
    rho_E = rho_E_alle{n_plot,i_plot};
    plot_dichtematrix(rho_E,Zellen,L,i_plot,n_plot);
end

%% --------- LOKALE FUNKTIONEN (am Dateiende belassen) ---------
function [Zellen,idmap] = rauten_gitter(K,N_xi)
    t_kanten = linspace(0,1,K+1);
    z_kanten = linspace(-1,1,N_xi+1);
    Raute = polyshape([0,.5,1,.5],[0,1,0,-1]);
    Zellen = struct('k',{},'j',{},'v',{},'M',{},'S',{});
    idmap = zeros(K,N_xi);

    for k=1:K
        for j=1:N_xi
            tL=t_kanten(k); tR=t_kanten(k+1);
            zL=z_kanten(j); zR=z_kanten(j+1);
            Rechteck=polyshape([tL,tR,tR,tL],[zL,zL,zR,zR]);
            Schnitt=intersect(Rechteck,Raute);
            if area(Schnitt)<1e-12
                continue
            end
            [t,z]=boundary(Schnitt);
            v=[t(:),z(:)];
            v=v(all(isfinite(v),2),:);
            if norm(v(end,:)-v(1,:))<1e-12
                v(end,:)=[];
            end
            if size(v,1)<3
                continue
            end
            nxt=[2:size(v,1),1];
            orientierung=sum(v(:,1).*v(nxt,2)-v(nxt,1).*v(:,2));
            if orientierung<0
                v=flipud(v); % Gegen den Uhrzeigersinn: richtige aeussere Normale
            end
            [Mlok,Slok]=zellen_integrale(v,tL,tR);
            idx=numel(Zellen)+1;
            Zellen(idx).k=k;
            Zellen(idx).j=j;
            Zellen(idx).v=v;
            Zellen(idx).M=Mlok;
            Zellen(idx).S=Slok;
            idmap(k,j)=idx;
        end
    end
end

function [Mlok,Slok]=zellen_integrale(v,tL,tR)
    % Exakte Integration fuer quadratische Produkte linearer DG-Basisfunktionen
    % ueber eine konvexe Polygonzelle mittels Dreieckszerlegung.
    Mlok=zeros(2,2);
    Slok=zeros(2,2);
    grad=[-1;1]/(tR-tL);
    bary=[2/3,1/6,1/6;1/6,2/3,1/6;1/6,1/6,2/3];
    for q=2:size(v,1)-1
        tri=v([1,q,q+1],:);
        d1=tri(2,:)-tri(1,:);
        d2=tri(3,:)-tri(1,:);
        Flaeche=abs(d1(1)*d2(2)-d1(2)*d2(1))/2;
        for a=1:3
            punkt=bary(a,:)*tri;
            ell=[(tR-punkt(1))/(tR-tL);(punkt(1)-tL)/(tR-tL)];
            Mlok=Mlok+(Flaeche/3)*(ell*ell.');
            Slok=Slok+(Flaeche/3)*(grad*ell.');
        end
    end
end

function T=rand_abbildung(R,Pin,Pout)
    % Gesucht: eingehende Anteile aus R*(u_out+u_in)=0.
    Aein=R*Pin; % 8x16, Rang 8 fuer ideale Andreev-Bedingungen
    if rank(Aein,1e-9)<8
        error('Randbedingung bestimmt die einlaufenden Charakteristiken nicht eindeutig.');
    end
    T=Pout-pinv(Aein)*(R*Pout);
end

function [K1,K2]=montiere_cutcell(Zellen,idmap,C1,C1p,C1m,C2,C2p,C2m,V1,V2,Rand)
    nCells=numel(Zellen);
    Nu=32*nCells;
    K1=spalloc(Nu,Nu,1000*Nu);
    K2=spalloc(Nu,Nu,1000*Nu);
    Kmax=size(idmap,1);
    Jmax=size(idmap,2);
    tgrid=linspace(0,1,Kmax+1);
    zgrid=linspace(-1,1,Jmax+1);
    tol=1e-8;
    gau=[(1-1/sqrt(3))/2,(1+1/sqrt(3))/2];

    for z=1:nCells
        k=Zellen(z).k; j=Zellen(z).j;
        tL=tgrid(k); tR=tgrid(k+1);
        M=Zellen(z).M; S=Zellen(z).S;
        vs=Zellen(z).v;

        %% Volumenanteile von Gleichung (10) und Gleichung (11)
        for ii=1:2
            zeilen=(z-1)*32+(ii-1)*16+(1:16);
            for mm=1:2
                spalten=(z-1)*32+(mm-1)*16+(1:16);
                K1(zeilen,spalten)=K1(zeilen,spalten) ...
                    + 1i/2*S(ii,mm)*C1 + M(ii,mm)*V1;
                K2(zeilen,spalten)=K2(zeilen,spalten) ...
                    - 1i/2*S(ii,mm)*C2 + M(ii,mm)*V2;
            end
        end

        %% Alle Polygonkanten: innere DG/FV-Grenzen oder Andreev-Rand
        for e=1:size(vs,1)
            e2=mod(e,size(vs,1))+1;
            p0=vs(e,:); p1=vs(e2,:);
            dp=p1-p0;
            laenge=norm(dp);
            if laenge<1e-12, continue; end
            normal=[dp(2),-dp(1)]/laenge;
            a1=normal(1)/2+normal(2);
            a2=normal(1)/2-normal(2);

            % Upwind-Splits der normalen Transportoperatoren
            G1p=max(a1,0)*C1p+min(a1,0)*C1m;
            G1m=max(a1,0)*C1m+min(a1,0)*C1p;
            G2p=max(a2,0)*C2p+min(a2,0)*C2m;
            G2m=max(a2,0)*C2m+min(a2,0)*C2p;

            % Nachbar durch gemeinsame vertikale/horizontale Rasterkante
            kn=0; jn=0;
            if abs(dp(1))<tol
                if abs(p0(1)-tL)<tol, kn=k-1; jn=j; end
                if abs(p0(1)-tR)<tol, kn=k+1; jn=j; end
            elseif abs(dp(2))<tol
                if abs(p0(2)-zgrid(j))<tol, kn=k; jn=j-1; end
                if abs(p0(2)-zgrid(j+1))<tol, kn=k; jn=j+1; end
            end
            nb=0;
            if kn>=1 && kn<=Kmax && jn>=1 && jn<=Jmax
                nb=idmap(kn,jn);
            end

            if nb==0
                % Physikalische schräge Rautenkante klassifizieren
                pm=(p0+p1)/2;
                x1=pm(1)+pm(2)/2;
                x2=pm(1)-pm(2)/2;
                T1=eye(16); T2=eye(16);
                if abs(x1)<tol
                    T1=Rand.T1L;
                elseif abs(x1-1)<tol
                    T1=Rand.T1R;
                elseif abs(x2)<tol
                    T2=Rand.T2L;
                elseif abs(x2-1)<tol
                    T2=Rand.T2R;
                else
                    error('Kante nicht als Nachbar oder physikalischer Rand identifiziert.');
                end
            end

            for q=1:2
                p=p0+gau(q)*dp;
                gewicht=laenge/2;
                ell=[(tR-p(1))/(tR-tL);(p(1)-tL)/(tR-tL)];
                if nb>0
                    tLn=tgrid(kn); tRn=tgrid(kn+1);
                    ellN=[(tRn-p(1))/(tRn-tLn);(p(1)-tLn)/(tRn-tLn)];
                end
                for ii=1:2
                    zeilen=(z-1)*32+(ii-1)*16+(1:16);
                    for mm=1:2
                        spalten=(z-1)*32+(mm-1)*16+(1:16);
                        if nb>0
                            K1(zeilen,spalten)=K1(zeilen,spalten) ...
                                -1i*gewicht*ell(ii)*ell(mm)*G1p;
                            K2(zeilen,spalten)=K2(zeilen,spalten) ...
                                +1i*gewicht*ell(ii)*ell(mm)*G2p;
                            spN=(nb-1)*32+(mm-1)*16+(1:16);
                            K1(zeilen,spN)=K1(zeilen,spN) ...
                                -1i*gewicht*ell(ii)*ellN(mm)*G1m;
                            K2(zeilen,spN)=K2(zeilen,spN) ...
                                +1i*gewicht*ell(ii)*ellN(mm)*G2m;
                        else
                            % Kein Nachbar: Charakteristischer Andreev-Randflux
                            K1(zeilen,spalten)=K1(zeilen,spalten) ...
                                -1i*gewicht*ell(ii)*ell(mm)*(G1p+G1m*T1);
                            K2(zeilen,spalten)=K2(zeilen,spalten) ...
                                +1i*gewicht*ell(ii)*ell(mm)*(G2p+G2m*T2);
                        end
                    end
                end
            end
        end
    end
end

function u=normiere_rho(u,Zellen,idmap,K,N_xi,L)
    % Normierung ueber die Diagonale x1=x2, d.h. z=0.
    twerte=linspace(0,1,201);
    werte=zeros(size(twerte));
    j0=(N_xi+1)/2;
    for p=1:numel(twerte)
        t=twerte(p);
        k=min(K,floor(t*K)+1);
        z=idmap(k,j0);
        if z==0, error('Kein aktiver Zustand bei xi=0.'); end
        tL=(k-1)/K; tR=k/K;
        ell=[(tR-t)/(tR-tL);(t-tL)/(tR-tL)];
        v=ell(1)*u((z-1)*32+(1:16)) ...
         +ell(2)*u((z-1)*32+(17:32));
        werte(p)=trace(reshape(v,4,4));
    end
    normwert=L*trapz(twerte,werte);
    if abs(normwert)<1e-14*L
        warning('Dichte-Normierung praktisch null; Zustand ggf. unphysikalisch.');
        return
    end
    u=u/normwert;
end

function plot_dichtematrix(rho_E,Zellen,L,i,n)
    % Zellenmitten: Darstellung als gefaerbte Punkte im (X,xi)-Gebiet.
    C=numel(Zellen);
    X=zeros(C,1); xi=zeros(C,1);
    werte=complex(zeros(4,4,C));
    Kmax=max([Zellen.k]);
    tgrid=linspace(0,1,Kmax+1);
    for c=1:C
        v=Zellen(c).v;
        p=polyshape(v(:,1),v(:,2));
        [t,z]=centroid(p);
        X(c)=L*t*1e9;
        xi(c)=L*z*1e9;
        k=Zellen(c).k;
        ell1=(tgrid(k+1)-t)/(tgrid(k+1)-tgrid(k));
        ell2=1-ell1;
        werte(:,:,c)=ell1*rho_E(:,:,1,c)+ell2*rho_E(:,:,2,c);
    end
    figure('Name',sprintf('rho_E, Phasenindex %d, ABS-Zweig %d',i,n));
    for a=1:4
        for b=1:4
            subplot(4,4,(a-1)*4+b);
            scatter(X,xi,24,squeeze(real(werte(a,b,:))),'filled');
            axis tight; axis equal;
            title(sprintf('Re rho_{%d%d}',a,b));
            if a==4, xlabel('X [nm]'); end
            if b==1, ylabel('xi [nm]'); end
        end
    end
    sgtitle(sprintf('DG(X)+FV(xi): rho_E, Phase %d, ABS-Zweig %d',i,n));
end
