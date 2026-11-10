function driftoperator = Potential(Ort, p)
driftoperator = zeros(size(Ort));

for i = 1 : length(Ort(1,:))
    for j = 1 : length(Ort(:,1))
        if (Ort(j, i) <= p.L_chi/2-p.bar_w/2 && Ort(j, i) >= p.L_chi/2-p.bar_g-p.bar_w/2 || Ort(j, i) >= p.L_chi/2+p.bar_w/2 && Ort(j, i) <= p.L_chi/2+p.bar_g+p.bar_w/2)
            driftoperator(j, i) = max(p.V_pot);
        else
            driftoperator(j, i) = min(p.V_pot);
        end
    end
end

end