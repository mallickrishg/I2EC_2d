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
mu = 1;% in MPa

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
BC([iright;iright]) = 0;
% bottom
BC([ibot;ibot]) = 0;
% left top half
BC([ileft1;ileft1]) = [boundary.Vpl(ileft1);boundary.Vpl(ileft1).*0];
% left bot half
BC([ileft2;ileft2]) = 0;

%% compute displacement & traction kernels for boundary
[Kdd,Kdn,Knd,Knn] = geometry.computeFullTractionKernels(boundary,boundary);
[Gdx,Gdz,Gnx,Gnz] = computeDisplacementKernels(boundary,boundary.xc);
