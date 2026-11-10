function [] = plotW(WW,n, nm)
if nm==2
    figure(n)
    tiledlayout(4,2)
    hold off; 
    nexttile
    imagesc(full(squeeze(real(WW(1,1,:,:)))));
    colorbar
    title('Real 11')

    nexttile
    imagesc(full(squeeze(real(WW(1,2,:,:)))));
    colorbar
    title('Real 12')

    nexttile
    imagesc(full(squeeze(real(WW(2,1,:,:)))));
    colorbar
    title('Real 21')

    nexttile
    imagesc(full(squeeze(real(WW(2,2,:,:)))));
    colorbar
    title('Real 22')

    nexttile
    imagesc(full(squeeze(imag(WW(1,1,:,:)))));
    colorbar
    title('Imag 11')

    nexttile
    imagesc(full(squeeze(imag(WW(1,2,:,:)))));
    colorbar
    title('Imag 12')

    nexttile
    imagesc(full(squeeze(imag(WW(2,1,:,:)))));
    colorbar
    title('Imag 21')

    nexttile
    imagesc(full(squeeze(imag(WW(2,2,:,:)))));
    colorbar
    title('Imag 22')

    drawnow;
end
end