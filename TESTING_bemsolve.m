% Script to setup a BEM problem that will later be used to construct
% displacement & stress kernels for a given mesh.
% We will apply displacement and traction BCs on different sections of the
% mesh. Need to figure out how to parameterize this for the user.
% 
% Author:
% Rishav Mallick, JPL, 2023

clear
addpath functions/
import geometry.*

% Elastic parameters (homogenous medium)
nu = 0.25;% Poisson's ratio
mu = 1e3;% in MPa

earthModel = geometry.LDhs(mu,nu);
rcv = geometry.receiver('inputs/testing2d.seg',earthModel);
boundary = geometry.receiver('inputs/boundary2d.seg',earthModel);

figure(1),clf
plotpatch2d(boundary,boundary.xc(:,2)), hold on
plotpatch2d(rcv,rcv.xc(:,2))
axis tight equal

%% compute source displacements & tractions at the boundary
[Kdd,Kdn,~,~] = geometry.computeFullTractionKernels(rcv,boundary);
[Gdx,Gdz,~,~] = geometry.computeDisplacementKernels(rcv,boundary.xc);

source = 1;
% calculate source tractions
td = Kdd*source;
tn = Kdn*source;
% calculate source displacements
ux = Gdx*source;
uz = Gdz*source;

% boundary condition vector
BC = zeros(2*boundary.N,1);
% boundary conditions we will use are as follows:
% right side: ux,uz = 0
% bottom: td,tn = 0
% left top half: ux = 1, uz = 0
% left bot half: td,tn = 0
iright = boundary.xc(:,1) == 10e3;
ibot = boundary.xc(:,2) == -10e3;
ileft1 = boundary.xc(:,1) == -10e3 & boundary.xc(:,2) >= -4e3;
ileft2 = boundary.xc(:,1) == -10e3 & boundary.xc(:,2) < -4e3;

% right side
BC([iright;iright]) = 0 - [ux(iright);uz(iright)];
% bottom
BC([ibot;ibot]) = 0 - [td(ibot);tn(ibot)];
% left top half
BC([ileft1;ileft1]) = [boundary.Vpl(ileft1)-ux(ileft1); 0-uz(ileft1)];
% left bot half
BC([ileft2;ileft2]) = 0 - [td(ileft2);tn(ileft2)];

%% compute displacement & traction kernels for boundary
[Kdd,Kdn,Knd,Knn] = geometry.computeFullTractionKernels(boundary,boundary);
[Gdx,Gdz,Gnx,Gnz] = geometry.computeDisplacementKernels(boundary,boundary.xc);

Ktraction = [Kdd,Knd;Kdn,Knn];
Kdisp = [Gdx,Gnx;Gdz,Gnz];

kernel = zeros(boundary.N*2,boundary.N*2);

kernel([iright;iright],:) = Kdisp([iright;iright],:);
kernel([ibot;ibot],:) = Ktraction([ibot;ibot],:);
kernel([ileft1;ileft1],:) = Kdisp([ileft1;ileft1],:);
kernel([ileft2;ileft2],:) = Ktraction([ileft2;ileft2],:);

%% solve boundary value problem

bem_sol = kernel\BC;

%% plot internal displacements
nx = 100;
nz = nx/2;

x = linspace(-9.99,9.99,nx).*1e3;
z = linspace(-9.99,0,nz).*1e3;
[X,Z] = meshgrid(x,z);
obs_plot = [X(:), Z(:)];

[Gdx,Gdz,Gnx,Gnz] = geometry.computeDisplacementKernels(boundary,obs_plot);

ux = [Gdx,Gnx]*bem_sol;
uz = [Gdz,Gnz]*bem_sol;

figure(10),clf
nskip = 8;
toplot = sqrt(ux.^2 + uz.^2);
pcolor(x,z,reshape(toplot,nz,nx)), shading interp, hold on
quiver(X(1:nskip:end),Z(1:nskip:end),ux(1:nskip:end)',uz(1:nskip:end)','w-','LineWidth',1)
axis tight equal
cb=colorbar;cb.Label.String = 'Displacement';
clim([0 1])
colormap(ttscm('berlin',20))
set(gca,'FontSize',20,'LineWidth',1.5,'TickDir','both')





















