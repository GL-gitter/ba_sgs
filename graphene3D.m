close all
clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Graphene 3D Dispersion Relation     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_real = 2.46;
t = -3.033;
s=0.129;
x=linspace(-3*pi/(sqrt(3)*a_real) ,3*pi/(sqrt(3)*a_real),1000);
y=linspace(-3.5*pi/(2*a_real),5*pi/(3*a_real),1000);
[a,b]=meshgrid(x,y);
E_bonding= t*sqrt(1+4*cos(a*a_real*sqrt(3)/2).*cos(b*a_real/2)...
    +4*(cos(b*a_real/2)).^2)./(1+s*sqrt(1+4*...
    cos(a*a_real*sqrt(3)/2).*cos(b*a_real/2)...
    +4*(cos(b*a_real/2)).^2));
E_antibonding=-t*sqrt(1+4*cos(a*a_real*sqrt(3)/2)...
    .*cos(b*a_real/2)+4*(cos(b*a_real/2)).^2)./...
    (1-s*sqrt(1+4*cos(a*a_real*sqrt(3)/2)...
    .*cos(b*a_real/2)+4*(cos(b*a_real/2)).^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Figure Setting            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
mesh(a,b,E_bonding)
xlabel('$k_{x}$','FontSize',20,'FontWeight','bold','Interpreter','latex')
ylabel('$k_{y}$','FontSize',20,'FontWeight','bold','Interpreter','latex')
zlabel('Energie (eV)','FontSize',20,'FontWeight','bold','Interpreter','latex')
hold on
mesh(a,b,E_antibonding)
c = colorbar;
c.Color = 'k';


set(gca, 'Color', 'w');   % Setzt den Hintergrund der Zeichenfläche (Achsenbereich) auf Weiß
set(gcf, 'Color', 'w');   % Setzt den Hintergrund des Figure-Fensters (äußerer Rahmen) auf Weiß
set(gca, 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k'); % Achsen und Skalenstriche auf Schwarz
set(gca, 'GridColor', 'k'); % Falls Gitter vorhanden sind