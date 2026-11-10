function [] = plotf(A)
    he = size(squeeze(A));
    if ndims(squeeze(A))==2 & he(2)==1
        plot(real(squeeze(A)))
    elseif ndims(squeeze(A))==2
        imagesc(real(squeeze(A)))
    end
end

