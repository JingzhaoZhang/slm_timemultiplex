function [loss, df ] = SourceFunObj( phase, source, z, Nx, Ny, thresholdh, thresholdl, maskFun, fresnelKernelFun, useGPU)
%FUNOBJ Summary of this function goes here
%   Detailed explanation goes here



loss = 0;

[nx, ~, ~] = size(phase);
numplex = nx/Nx;

df = zeros(nx, Ny);


imageInten = zeros(Nx, Ny);
imagez = zeros(nx, Ny);
temph = zeros(Nx, Ny);
templ = zeros(Nx, Ny);
if useGPU
    imageInten = gpuArray(imageInten);
    imagez =  gpuArray(imagez);
    temph = gpuArray(temph);
    templ = gpuArray(templ);
end

for i = 1 : numel(z)
    imageInten(:)= 0;
    imagez(:)= 0;
    temph(:)= 0;
    templ(:)= 0;
    
    HStack = fresnelKernelFun(i,i);
    mask = maskFun(i,i);    
    
    for n = 1:numplex
        objectField = gpuArray(source.*exp(1i * phase(Nx*(n-1)+ 1:Nx*n, :)));
        imagez(Nx*(n-1)+ 1:Nx*n, :) = fftshift(fft2(objectField .* HStack));
        imageInten = imageInten + abs(imagez(Nx*(n-1)+ 1:Nx*n, :).^2);
    end
    
    maskh = mask .* (imageInten < thresholdh);
    maskl = (1-mask) .* (imageInten > thresholdl);
    
    diffh = maskh .* (imageInten - thresholdh);
    diffl = maskl .* (imageInten - thresholdl);
    
    for n = 1:numplex
        temph = imagez(Nx*(n-1)+ 1:Nx*n, :).*diffh;
        temph = conj(HStack).*(Nx*Ny*ifft2(...
            ifftshift(temph)));
        templ = imagez(Nx*(n-1)+ 1:Nx*n, :).*diffl;
        templ = conj(HStack).*(Nx*Ny*ifft2(...
            ifftshift(templ)));
        df(Nx*(n-1)+ 1:Nx*n, :) =  df(Nx*(n-1)+ 1:Nx*n, :) + gather(temph+templ);
    end
%     templ = Nx*Ny *abs(HStack.^2).* (objectField.*diffl);
%     temph = Nx*Ny *abs(HStack.^2).* (objectField.*diffh);

    loss = loss + sum(sum(diffh.^2 + diffl.^2));     
    %df = df +  temph + templ;
    %clear HStack mask imagez imageInten maskh maskl diffh diffl temph templ
end
%df = real(source.*(-real(df).*sin(phase)) + imag(df).*cos(phase));

for n = 1:numplex
dfphase(Nx*(n-1)+ 1:Nx*n, :) = source.*(- real(df(Nx*(n-1)+ 1:Nx*n, :)...
    ).*sin(phase(Nx*(n-1)+ 1:Nx*n, :))...
    + imag(df(Nx*(n-1)+ 1:Nx*n, :)) .* cos(phase(Nx*(n-1)+ 1:Nx*n, :)));
end


df = gather(real(dfphase));
% loss = real(loss);
% df(1:Nx*Ny) = real(dfphase(:)) * ratio1;
% df(Nx*Ny+1:end) = real(dfsource(:))* ratio2;
loss = gather(real(loss));
end

