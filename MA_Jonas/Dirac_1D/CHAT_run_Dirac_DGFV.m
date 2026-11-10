function CHAT_run_Dirac_DGFV()
% RUN_DIRAC_DGFV
% Schritt 0: freier 1D-Dirac-Kanal (ky=0, Delta=0, B=0).
% Numerik: DG vom Grad 1 in X und Finite Volumes in xi,
%          periodische Randbedingungen, RK4 in der Zeit.
%
% Ablage: Diese Datei neben dem Ordner DG_code speichern.
%         Es werden drei vorhandene DG-Core-Funktionen verwendet:
%         JacobiGL.m, Vandermonde1D.m, Dmatrix1D.m.
%
% Die periodische Rechteck-Testgeometrie in (X,xi) ist ein
% numerischer Testfall, noch keine Kontaktgeometrie des SGS-Kontakts.

    % (1) Vorhandene DG-Basisfunktionen zugreifbar machen.
    here = fileparts(mfilename('fullpath'));
    core = fullfile(here, 'DG_code', 'core');
    if ~isfolder(core)
        core = fullfile(here, 'core');
    end
    assert(isfolder(core), ...
        'Ordner DG_code/core nicht gefunden. Datei neben DG_code ablegen.');
    addpath(core);

    % (2) Physikalische und numerische Parameter.
    vF = 1.0e6;          % Fermi-Geschwindigkeit [m/s]
    L = 1.0e-6;          % Laenge des periodischen Testbereichs [m]
    K = 24;              % DG-Elemente in X
    Np = 2;              % 2 DG-Knoten je Element -> Polynomgrad 1
    Nxi = 31;            % FV-Zellen in xi (ungerade -> xi=0 ist Zellmitte)
    hX = L/K;
    hXi = L/Nxi;

    % (3) DG-Matrizen: diese Routinen stammen aus DG_code/core.
    r = JacobiGL(0, 0, Np-1);
    V = Vandermonde1D(Np-1, r);
    Dr = Dmatrix1D(Np-1, r, V);  % Ableitung bzgl. Referenz-x
    invV = inv(V);
    M = invV.' * invV;            % DG-Massenmatrix
    invM = M \ eye(Np);

    % (4) Dirac-Transportmatrizen fuer u=vec(rho).
    sigmaX = [0,1;1,0];
    A = vF * sigmaX;
    I2 = eye(2);
    AX = 0.5 * (kron(I2,A) + kron(A.',I2));
    Axi = kron(I2,A) - kron(A.',I2);

    % Matrix-Upwind: positive/negative charakteristische Geschwindigkeiten.
    [AXplus, AXminus] = split_flux(AX);
    [Axiplus, Aximinus] = split_flux(Axi);

    % (5) Raeumliche Koordinaten und Anfangszustand.
    xi = ((1:Nxi) - (Nxi+1)/2) * hXi;
    Xnodes = zeros(Np,K);
    for k = 1:K
        Xnodes(:,k) = (k-1)*hX + (r+1)*hX/2;
    end

    x0 = 0.27*L;           % Mittelpunkt der Anfangswelle
    sigma = 0.07*L;        % Wellenpaketbreite
    sPlus = [1;1]/sqrt(2); % Rechtslaeufiger Dirac-Pseudospin
    P = sPlus*sPlus';       % 2x2 Dichtematrix des Pseudospins
    pvec = P(:);           % vier Komponenten: vec(rho)

    % U(a,i,k,j) = Komponente a am DG-Knoten i, im Element k,
    %                und in der FV-Zelle j von xi.
    U = zeros(4,Np,K,Nxi);
    for j = 1:Nxi
        for k = 1:K
            for i = 1:Np
                X = Xnodes(i,k);
                % rho_ana(X,xi,0) aus dem analytischen Wellenpaket:
                f = exp(-((X-x0)/sigma)^2 - (xi(j)/(2*sigma))^2);
                U(:,i,k,j) = f*pvec;
            end
        end
    end

    % (6) Zeitentwicklung der DISKRETISIERTEN Differentialgleichung.
    % Die analytische Loesung wird dabei NICHT verwendet.
    tEnd = 0.25*L/vF;
    dt = 0.15*min(hX/vF, hXi/(2*vF));
    Nt = ceil(tEnd/dt);
    dt = tEnd/Nt;

    mass0 = particle_number(U(:,:,:, (Nxi+1)/2), M, hX);
    for n = 1:Nt
        k1 = rhs_dirac(U, AX, AXplus, AXminus, Axi, Axiplus, Aximinus, Dr, invM, hX, hXi);
        k2 = rhs_dirac(U+dt*k1/2, AX, AXplus, AXminus, Axi, Axiplus, Aximinus, Dr, invM, hX, hXi);
        k3 = rhs_dirac(U+dt*k2/2, AX, AXplus, AXminus, Axi, Axiplus, Aximinus, Dr, invM, hX, hXi);
        k4 = rhs_dirac(U+dt*k3, AX, AXplus, AXminus, Axi, Axiplus, Aximinus, Dr, invM, hX, hXi);
        U = U + dt*(k1+2*k2+2*k3+k4)/6;
    end

    % (7) Vergleich der numerischen Dichte bei xi=0 mit der Analytik.
    j0 = (Nxi+1)/2;
    xPlot = Xnodes(:)*1e6;  % Mikrometer
    nNum = zeros(Np,K);
    for k = 1:K
        for i = 1:Np
            R = reshape(U(:,i,k,j0),2,2);
            nNum(i,k) = real(trace(R));
        end
    end
    nAna = exp(-((Xnodes-x0-vF*tEnd)/sigma).^2);
    relErr = norm(nNum(:)-nAna(:))/norm(nAna(:));
    massEnd = particle_number(U(:,:,:,j0), M, hX);

    figure;
    plot(xPlot, nNum(:), 'o', 'DisplayName','DG(X)+FV(xi)');
    hold on;
    plot(xPlot, nAna(:), '-', 'LineWidth',1.5, 'DisplayName','Analytische Loesung');
    hold off; grid on;
    xlabel('X [\mum]');
    ylabel('n(X,t) = Tr(\rho(X,0,t))');
    title('Freier Dirac-Transport: erster numerischer Test');
    legend('Location','best');

    fprintf('Zeit tEnd = %.3f ps\n', tEnd*1e12);
    fprintf('Relativer Fehler der Dichte = %.4g\n', relErr);
    fprintf('Relative Aenderung der Gesamtteilchenzahl = %.4g\n', ...
        (massEnd-mass0)/mass0);
    fprintf(['Hinweis: Erstes Beispiel testet vor allem den rechtslaufenden ', ...
             'Schwerpunkttransport.\n']);
    fprintf(['Fuer eine vollstaendige Validierung muss spaeter auch ', ...
             'der xi-Transport angeregt werden.\n']);
end

function [Aplus, Aminus] = split_flux(A)
% Zerlegung A = Aplus + Aminus nach Eigenwert-Vorzeichen.
    [V,D] = eig(A);
    lambda = real(diag(D));
    Aplus = real(V*diag(max(lambda,0))*V');
    Aminus = real(V*diag(min(lambda,0))*V');
end

function dU = rhs_dirac(U, AX, AXplus, AXminus, Axi, Axiplus, Aximinus, Dr, invM, hX, hXi)
% Räumlicher DG(X)- und FV(xi)-Operator mit periodischen Flüssen.
    [~,Np,K,Nxi] = size(U);
    dU = zeros(size(U));

    % DG in X: Ableitung innerhalb eines Elements und Upwind-Grenzfluss.
    for j = 1:Nxi
        for k = 1:K
            km = mod(k-2,K)+1;
            kp = mod(k,K)+1;
            Uk = U(:,:,k,j);
            Fk = AX*Uk;
            Fleft = AXplus*U(:,Np,km,j) + AXminus*Uk(:,1);
            Fright = AXplus*Uk(:,Np) + AXminus*U(:,1,kp,j);

            corr = zeros(4,Np);
            corr(:,1) = Fleft - Fk(:,1);
            corr(:,Np) = corr(:,Np) + Fk(:,Np) - Fright;

            dU(:,:,k,j) = -(2/hX)*(Fk*Dr.') + (2/hX)*(corr*invM.');
        end
    end

    % FV in xi: numerischer Fluss über jede Zellgrenze.
    for j = 1:Nxi
        jm = mod(j-2,Nxi)+1;
        jp = mod(j,Nxi)+1;
        uj = reshape(U(:,:,:,j),4,[]);
        um = reshape(U(:,:,:,jm),4,[]);
        up = reshape(U(:,:,:,jp),4,[]);
        Fright = Axiplus*uj + Aximinus*up;
        Fleft = Axiplus*um + Aximinus*uj;
        dU(:,:,:,j) = dU(:,:,:,j) - reshape((Fright-Fleft)/hXi,4,Np,K);
    end
end

function N = particle_number(Ucenter, M, hX)
% Integral ueber X der Dichte bei xi=0 (nur Test fuer dieses Wellenpaket).
    [~,Np,K] = size(Ucenter);
    N = 0;
    w = M*ones(Np,1); % Gewicht fuer Integral über DG-Element
    for k = 1:K
        n = zeros(Np,1);
        for i = 1:Np
            R = reshape(Ucenter(:,i,k),2,2);
            n(i) = real(trace(R));
        end
        N = N + hX/2 * (w.'*n);
    end
end
