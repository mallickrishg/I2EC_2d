clear
dip = 20*pi/180;
nx = 1000;
ny = 500;
x = linspace(-2,4,nx)';
z = linspace(-4*tan(dip),-0.001,ny)';
[X,Z] = meshgrid(x,z);
[Vx,Vz] = compute_cornerflow_velocityfield(dip,X(:),Z(:));

figure(1),clf
index = Z(:)+X(:).*tan(dip)<=0;
nskip = 410;
toplotx = Vx(:,1);
toplotz = Vz(:,1);
toplotx(index) = nan;
toplotz(index) = nan;
pcolor(x,z,reshape(sqrt(toplotx.^2 + toplotz.^2),ny,nx)), shading interp, hold on
quiver(X(1:nskip:end)',Z(1:nskip:end)',toplotx(1:nskip:end),toplotz(1:nskip:end),'k-','Linewidth',1)

index = Z(:)+X(:).*tan(dip)>=0;
toplotx = Vx(:,2);
toplotz = Vz(:,2);
toplotx(index) = nan;
toplotz(index) = nan;
pcolor(x,z,reshape(sqrt(toplotx.^2 + toplotz.^2),ny,nx)), shading interp, hold on
quiver(X(1:nskip:end)',Z(1:nskip:end)',toplotx(1:nskip:end),toplotz(1:nskip:end),'w-','Linewidth',1)
axis tight equal
colormap turbo(10)
colorbar
clim([0 1])