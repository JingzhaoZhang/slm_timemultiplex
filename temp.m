%%
% f = figure();
% n = numel(z);
% for i = 1:n
%     i
%     imagesc(Ividmeas(:,:,i));colormap gray;axis off;colorbar;
%     caxis([0,1]);
%     %colorbar;
%  
% end
% close(f);

%%
close all
filename = 'PointCloud_819_800by600_points_gaussian_stack';
load(['data/' filename '.mat']);
figure();
for i = 1:81
    i;
    imagesc(Ividmeas(:,:,i));colormap gray;
end

figure()
slopeX = (388- 553)/(67-1);
slopeY = (426 - 189)/(67-1);
Ividmeas1 = zeros(size(Ividmeas));
for i = 1:81
    i;
    Im = circshift(Ividmeas(:,:,i), [ -floor((i-37)*slopeY), -floor((i-37)*slopeX)]);
    imagesc(Im);colormap gray;caxis([0,1]);
    pause(0.05);
    Ividmeas1(:,:,i) = Im;
end


%%
filename = 'optimization_819_600by800_points_gaussian_stack';
load(['data/' filename '.mat']);

figure();
for i = 1:81
    i;
    imagesc(Ividmeas(:,:,i));colormap gray;
end

figure()
slopeX = (429-592)/(68-1);
slopeY = (423-183)/(68-1);
Ividmeas2 = zeros(size(Ividmeas));

for i = 1:81
    i;
    Im = circshift(Ividmeas(:,:,i), [ -floor((i-39)*slopeY), -floor((i-39)*slopeX)]);
    imagesc(Im);colormap gray;caxis([0,1]);
    pause(0.05);
    Ividmeas2(:,:,i) = Im;
end




%%
x1 = 315; 
y1 = 594; 
x2 = 318;
y2 = 629;
t = 30;

v1 = squeeze(sum(sum(Ividmeas1(x1-t:x1+t, y1-t:y1+t, :))));
v1 = v1/max(v1);
v2 = squeeze(sum(sum(Ividmeas2(x2-t:x2+t, y2-t:y2+t, :))));
v2 = v2/max(v2);
figure();plot(1:numel(z), v1, 1:numel(z), v2);legend('Point Cloud', 'Optomization');

%%
x1 = 400; 
y1 = 300; 
x2 = 400;
y2 = 300;
t = 5;

v1 = squeeze(sum(sum(Ividmeas1(x1-t:x1+t, y1-t:y1+t, :))));
v1 = v1/max(v1);
v2 = squeeze(sum(sum(Ividmeas2(x2-t:x2+t, y2-t:y2+t, :))));
v2 = v2/max(v2);
figure();plot(1:numel(z), v1, 1:numel(z), v2);legend('Optomization', 'Point Cloud');

