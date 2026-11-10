function V = get_V(Ort, p)
V = zeros(size(Ort));

for j = 1 : length(Ort(:))
    if (Ort(j) <= p.L_chi/2-p.bar_w/2 && Ort(j) >= p.L_chi/2-p.bar_g-p.bar_w/2 || Ort(j) >= p.L_chi/2+p.bar_w/2 && Ort(j) <= p.L_chi/2+p.bar_g+p.bar_w/2)
        V(j) = 1.60343;
    else
        V(j) = 1.424;
    end
end

end
