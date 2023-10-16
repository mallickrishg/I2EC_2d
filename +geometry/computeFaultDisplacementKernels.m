function [Gdx,Gdz,Gnx,Gnz] = computeFaultDisplacementKernels(rcv,ox)

delta = 1e-9.*mean(ox(:,1));
xobs = ox(:,1);
zobs = ox(:,2);

Gdx = zeros(length(xobs),rcv.N);
Gdz = zeros(length(xobs),rcv.N);
Gnx = zeros(length(xobs),rcv.N);
Gnz = zeros(length(xobs),rcv.N);

for k = 1:rcv.N
    xobs_mod = xobs;
    zobs_mod = zobs;
    r = sqrt((xobs-rcv.xc(k,1)).^2 + (zobs-rcv.xc(k,2)).^2);
    xobs_mod(r==0) = xobs_mod(r==0) + rcv.nv(r==0,1)*delta;
    zobs_mod(r==0) = zobs_mod(r==0) + rcv.nv(r==0,2)*delta;
    
    [Disp] = geometry.LDdispHS(xobs_mod,zobs_mod,rcv.xc(k,1),rcv.xc(k,2),rcv.W(k)/2,...
        -deg2rad(rcv.dip(k)),1,0,rcv.earthModel.nu);
    
    Gdx(:,k) = Disp(:,1);
    Gdz(:,k) = Disp(:,2);
    
    [Disp] = geometry.LDdispHS(xobs_mod,zobs_mod,rcv.xc(k,1),rcv.xc(k,2),rcv.W(k)/2,...
        -deg2rad(rcv.dip(k)),0,1,rcv.earthModel.nu);
    
    Gnx(:,k) = Disp(:,1);
    Gnz(:,k) = Disp(:,2);
end


end