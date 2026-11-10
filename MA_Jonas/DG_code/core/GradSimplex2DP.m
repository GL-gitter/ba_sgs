function [dmodedr, dmodeds] = GradSimplex2DP(a,b,id,jd)
% function [dmodedr, dmodeds] = GradSimplex2DP(a,b,id,jd)
% Purpose: Return the derivatives of the modal basis (id,jd)
% on the 2D simplex at (a,b).

fa = JacobiP(a, 0, 0, id);
dfa = GradJacobiP(a, 0, 0, id);
gb = JacobiP(b, 2*id+1,0, jd); dgb = GradJacobiP(b, 2*id+1,0, jd);
dmodedr = dfa.*gb;

if(id>0)
    dmodedr = dmodedr.*((0.5*(1-b)).^(id-1));
end

dmodeds = dfa.*(gb.*(0.5*(1+a)));

if(id>0)
    dmodeds = dmodeds.*((0.5*(1-b)).^(id-1));
end

tmp = dgb.*((0.5*(1-b)).^id);

if(id>0)
    tmp = tmp-0.5*id*gb.*((0.5*(1-b)).^(id-1));
end

dmodeds = dmodeds+fa.*tmp;
% Normalize
dmodedr = 2^(id+0.5)*dmodedr; dmodeds = 2^(id+0.5)*dmodeds;
return;