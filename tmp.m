close all;
Ividmeas = Ividmeas/max(max(max(Ividmeas)));
n = size(Ividmeas, 3);
figure()
for i = 1:n
    imagesc(Ividmeas(:,:,i));colormap gray; colorbar; axis off;caxis([0,1])
end