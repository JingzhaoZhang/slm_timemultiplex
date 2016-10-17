resolutionScale = 1; % The demagnification scale of tubelens and objective. f_tube/f_objective
lambda = 0.532e-6;  % Wavelength
focal_SLM = 0.25; % focal length of the lens after slm.
psSLM = 8e-6;      % Pixel Size (resolution) at the scattered 3D region
Nx = 1920;       % Number of pixels in X direction
Ny = 1080;       % Number of pixels in Y direction


psXHolograph = lambda * focal_SLM/ psSLM / resolutionScale / Nx;      % Pixel Size (resolution) at the scattered 3D region
psYHolograph = lambda * focal_SLM/ psSLM / resolutionScale / Ny;      % Pixel Size (resolution) at the scattered 3D region

useGPU = 1;     % Use GPU to accelerate computation. Effective when Nx, Ny is large (e.g. 600*800).

std = Nx/8;
backAperture = 0.4;
z = [-100 : 2: 100] * 0.5e-4 ;   % Depth level requested in 3D region.
nfocus = floor(numel(z)/2) + 1;                % z(nfocus) denotes the depth of the focal plane.
thresholdh = 3e8;          % Intensity required to activate neuron.
thresholdl = 0;             % Intensity required to not to activate neuron.

%% Compute NA
cx=floor(Nx/2)+1;cy=floor(Ny/2)+1;
[us, vs]=ndgrid([1:Nx]-cx,[1:Ny]-cy);

masks = zeros(Nx, Ny, numel(z));
edge = 25;
r = z(end-edge) - z(edge);
dis = (us * psXHolograph).^2 + (vs * psYHolograph).^2;
figure()
for i = edge:numel(z)-edge
    masks(:,:,i) = 1 * (((z(i)- z(edge))^2 + dis) < r^2);
    imagesc(masks(:,:,i));title(int2str(i))
    pause(0.1)
end

save('half_sphere_mask.mat', 'masks')