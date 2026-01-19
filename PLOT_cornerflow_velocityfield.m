clear
addpath functions/

% specify subduction zone dip angle
dip = 25*pi/180;

% mesh points to plot
nx = 1000;
ny = 500;
x = linspace(-2,4,nx)';
z = linspace(-4*tan(dip),-0.001,ny)';
[X,Z] = meshgrid(x,z);

% compute velocity field in continental & oceanic mantle
% Vx,Vz are both [Nx2] matrices where
% Vx(:,1) - continental arc solution
% Vx(:,2) - oceanic mantle solution
% similarly for Vz

n_powerlaw = 3;
% [Vx,Vz] = compute_cornerflow_velocityfield(dip,X(:),Z(:));
% function that computes velocity field for n=1 or n=3 rheologies 
[Vx,Vz] = computeCornerflowVelocityfield(dip,X(:),Z(:),n_powerlaw);

% Plot velocity field
figure(1),clf

% solution only applies to continental arc
index = Z(:)+X(:).*tan(dip)<=0;

nskip = 410;
toplotx = Vx(:,1);
toplotz = Vz(:,1);
toplotx(index) = nan;
toplotz(index) = nan;
pcolor(x,z,reshape(sqrt(toplotx.^2 + toplotz.^2),ny,nx)), shading interp, hold on
quiver(X(1:nskip:end)',Z(1:nskip:end)',toplotx(1:nskip:end),toplotz(1:nskip:end),'k-','Linewidth',1)

% solution only applies to oceanic arc
index = Z(:)+X(:).*tan(dip)>=0;

toplotx = Vx(:,2);
toplotz = Vz(:,2);
toplotx(index) = nan;
toplotz(index) = nan;
pcolor(x,z,reshape(sqrt(toplotx.^2 + toplotz.^2),ny,nx)), shading interp, hold on
quiver(X(1:nskip:end)',Z(1:nskip:end)',toplotx(1:nskip:end),toplotz(1:nskip:end),'w-','Linewidth',1)
axis tight equal
colormap turbo(10)
cb=colorbar; cb.Label.String = 'velocity';
clim([0 1])
ylim([-1.5 0.1]), xlim([-2,4])
xlabel('x'), ylabel('z')