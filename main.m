close all;clear all;


%% Setup params
% All length value has unit meter in this file.
% The 3d region is behind lens after SLM.
tag = '1080p_1016_onebyone10multiplex';

resolutionScale = 1; % The demagnification scale of tubelens and objective. f_tube/f_objective
lambda = 0.532e-6;  % Wavelength
focal_SLM = 0.25; % focal length of the lens after slm.
psSLM = 8e-6;      % Pixel Size (resolution) at the scattered 3D region
Nx = 1920;       % Number of pixels in X direction
Ny = 1080;       % Number of pixels in Y direction
numplex = 10;


psXHolograph = lambda * focal_SLM/ psSLM / resolutionScale / Nx;      % Pixel Size (resolution) at the scattered 3D region
psYHolograph = lambda * focal_SLM/ psSLM / resolutionScale / Ny;      % Pixel Size (resolution) at the scattered 3D region

useGPU = 1;     % Use GPU to accelerate computation. Effective when Nx, Ny is large (e.g. 600*800).

std = Nx/8;
backAperture = 0.4;
z0 = [-100 : 0.5: 100] * 0.5e-4 ;   % Depth level requested in 3D region.
nfocus = floor(numel(z0)/2) + 1;                % z(nfocus) denotes the depth of the focal plane.
z0 = z0(1:end-1);
thresholdh = 5e7;          % Intensity required to activate neuron.
thresholdl = -2e7;             % Intensity required to not to activate neuron.

%% Compute NA
cx=floor(Nx/2)+1;cy=floor(Ny/2)+1;
[us, vs]=ndgrid([1:Nx]-cx,[1:Ny]-cy);

dis = (us * psXHolograph).^2 + (vs * psYHolograph).^2;
us=us * psSLM; vs=vs * psSLM;
%us = gpuArray(us); vs = gpuArray(vs);
r=(us.^2+vs.^2).^(1/2);
Pupil= r<=backAperture/2;


%% Kernel Function


% kernelfun = @(i) GenerateFresnelPropagationStack(Nx, Ny, z(i)-z(nfocus), lambda,psXHolograph,psYHolograph, focal_SLM, useGPU);


%% Pick Source Initialization method

% The starting point. reshape(x0(1:Nx*Ny), [Nx, Ny]) encodes the phase on
% SLM in rads. Normally initialized to zeros. reshape(x0(1+Nx*Ny:end), [Nx, Ny])
% encodes the source intensity. Need to be nonnegative.
subsetportion = 10; % divisible by 400
intensity = 1;
source = sqrt(intensity) * ones(Nx, Ny);
tag = [tag, '_uniform'];
Ividmeas = zeros(Nx, Ny, numel(z0)/subsetportion);
source = source/(max(max(source)));
masks = zeros(Nx, Ny, numel(z0));
r = 0.005;

%figure()
for i = 1:numel(z0)
    masks(:,:,i) = 1 * (((z0(i)- z0(100))^2 + dis) < r^2 * (z0(i) > z0(100)));
    %imagesc(masks(:,:,i));title(int2str(i))
    %pause(0.1)
end

HStacks = zeros(Nx, Ny, numel(z0));
for i = 1 : numel(z0)
    HStacks(:,:,i) = GenerateFresnelPropagationStack(Nx, Ny, z0(i)-z0(nfocus), lambda, psXHolograph,psYHolograph, 0);
end



%%
tic;

for N = 1:100
    
    if N > numplex && mod(N, 3) == 0
        save(['optimization_' tag '.mat'], 'phase', 'z');
    end    
    offset = mod(N, subsetportion)+1;
    z = z0(offset:subsetportion:end);
    imdepths = z;
    ims = gpuArray(masks(:,:,offset:subsetportion:end));
    maskfun = @(i1, i2) ims(:,:,i1:i2);
    substack = gpuArray(HStacks(:,:,offset:subsetportion:end));
    kernelfun = @(i1, i2) substack(:,:,i1:i2);
    
    if N > 1 && N <= numplex
        Ividmeas = zeros(Nx, Ny, numel(z), 'gpuArray');
        for i = 1:numel(z)
            HStack = substack(:,:,i);
            for n = 1:N-1
                Ividmeas(:,:,i) = Ividmeas(:,:,i) + fresnelProp(x0((n-1)*Nx+1:Nx*n, :), source, HStack);
            end
        end
        Ividmeas = Ividmeas/max(Ividmeas(:));
    end
    
    if N <= numplex
        hologram = zeros(Nx, Ny);
        for i = 1:numel(imdepths)
            HStack = substack(:,:,i);
            target = max(ims(:,:,i) - Ividmeas(:,:,i), 0);
            %    target = ims(:,:,i);
            target_phase = randn(Nx, Ny);
            imagez = sqrt(target) .* exp(1i * target_phase);
            im =  ifft2(ifftshift(imagez))./HStack;
            hologram = hologram + im;
        end
        x0((N-1)*Nx+1:N*Nx, :) = gather(angle(hologram));
        maxiter = 10;
    else
        maxiter = 3;
    end
    %% Optimization phase 1
    % Scale the gradient of phase by 1 and the gradient of source by 0.
    % This makes sure that only the phase is updated in each iteration.
    
    %tic;
    f = @(x)SourceFunObj( x, source, z, Nx, Ny, thresholdh * (1/2+min(N/2, numplex/2)),...
        thresholdl* (1/2+min(N/2, numplex/2)), maskfun, kernelfun, useGPU);
    
    matlab_options = optimoptions('fmincon','GradObj','on', 'display', 'iter', ...
        'algorithm','interior-point','Hessian','lbfgs', 'MaxFunEvals', 200, 'MaxIter', maxiter,...
        'TolX', 1e-20, 'TolFun', 1e-12);
    lb = -inf(Nx*Ny, 1);
    ub = inf(Nx*Ny, 1);
    nonlcon = [];
    phase = fmincon(f,x0,[],[],[],[],lb,ub,nonlcon,matlab_options);
    phase = reshape(phase, [Nx * min(N, numplex), Ny]);
    %phase = reshape(x0, [Nx * N, Ny]);
    x0 = mod(phase, 2*pi);
    
    toc;

end

%% plot
%Ividmeas = zeros(Nx, Ny, numel(z0));
usenoGPU = 0;
figure();
for i = 1:20:numel(z0)
    HStack = GenerateFresnelPropagationStack(Nx,Ny,z0(i) - z0(nfocus), lambda, psXHolograph,psYHolograph, usenoGPU);
    imagez = 0;
    for n = 1:numplex
        imagez = imagez + fresnelProp(phase((n-1)*Nx+1:Nx*n, :), source, HStack);
    end
    %h = fspecial('gaussian', [4,4], 1);
    %imagez = imfilter(imagez, h,'replicate');
    
    
    %Ividmeas(:,:,i) = imagez;
    imagesc(imagez);colormap gray;title(sprintf('Distance z %d', z0(i)));colorbar;
    caxis([0, 2e8]);
    pause(0.01);
end
%prctile(Ividmeas(:), 99.99)
%prctile(Ividmeas(:), 99.9)
%max(Ividmeas(:))
save(['optimization_' tag '.mat'], 'source', 'phase', 'z');
%save(['source_phase_result_' tag '.mat'], 'source1', 'phase1', 'source2', 'phase2', 'hologram');

%save(['simultaneous_result_' tag '.mat'],  'source2', 'phase2', 'hologram');

