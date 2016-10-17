function [ mask ] = generateComplexMask( zi, Nx, Ny, ims, imdepths)
% GENERATEMASK targets is a list of triple tuple [x,y,z; ...]. (0,0, z')
% is at the center. z' and vector z is within the same scale. All have
% unit meter.
% z is a sequence of z distances.
% ps is the pixel size.
% Returns the i th mask. i <= numel(z). The ball around targets has value
% 1.
mask = zeros(Nx, Ny);
for i = 1:numel(imdepths)
    if zi == imdepths(i)
        mask = ims(:,:,i);
        break
    end
end

end

