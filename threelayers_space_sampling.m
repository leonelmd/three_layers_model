function [kx, kz, ky] = threelayers_space_sampling(x,z)
%%%function [kx, kz, ky, kya] = threelayers_space_sampling(x,z)
%%%Samples space for multilayer planar model

if length(x)==1 || length(z)==1
    kx = 1;
    kz = 1;
    ky = 1;
    return
end
dx = x(2)-x(1);
dz = z(2)-z(1);
Lx = length(x);
Lz = length(z);

if mod(Lx,2)==0 %even
    dkx = 2*pi/(Lx*dx);
    gridx = ((-pi/dx):dkx:(pi/dx-dkx)); %see Mesin & Merletti, 2008, p.664
else
    dkx = 2*pi/((Lx-1)*dx);
    gridx = ((-pi/dx):dkx:(pi/dx)); %see Mesin & Merletti, 2008, p.664
end

if mod(Lz,2)==0 %even
    dkz = 2*pi/(Lz*dz);
    gridz = ((-pi/dz):dkz:(pi/dz-dkz)); 
else
    dkz = 2*pi/((Lz-1)*dz);
    gridz = ((-pi/dz):dkz:(pi/dz)); 
end

[kx, kz] = meshgrid(gridx,gridz);
ky = sqrt(kx.^2+kz.^2);