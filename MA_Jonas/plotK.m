function [] = plotK(K,n, nm)
figure(n)
tiledlayout(3,3)
hold off; nexttile
for IM = 1:nm
    plot(squeeze(K.K0(IM,IM,:))); hold on;
end
hold off; nexttile
for IM = 1:nm
    plot(squeeze(K.K1(IM,IM,:))); hold on;
end
hold off; nexttile
for IM = 1:nm
    plot(squeeze(K.K2(IM,IM,:))); hold on;
end
hold off; nexttile
for IMO = 1:nm
    for IMI=1:nm
        if IMO~=IMI
            plot(squeeze(K.K0(IMO,IMI,:))); hold on;
        end

    end
end
hold off; nexttile
for IMO = 1:nm
    for IMI=1:nm
        if IMO~=IMI
            plot(squeeze(K.K1(IMO,IMI,:))); hold on;
        end

    end
end
hold off; nexttile
for IMO = 1:nm
    for IMI=1:nm
        if IMO~=IMI
            plot(squeeze(K.K2(IMO,IMI,:))); hold on;
        end

    end
end
hold off; nexttile
imagesc(squeeze(max(abs(K.K0),[],3)));

hold off; nexttile
imagesc(squeeze(max(abs(K.K1),[],3)));

hold off; nexttile
imagesc(squeeze(max(abs(K.K2),[],3)));    

drawnow;
end