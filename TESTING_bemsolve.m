
clear
addpath functions/
import geometry.*

% Elastic parameters (homogenous medium)
nu=0.25;% Poisson's ratio
mu=1;% in MPa

earthModel = geometry.LDhs(mu,nu);
rcv = geometry.receiver('inputs/testing2d.seg',earthModel);
boundary = geometry.receiver('inputs/boundary2d.seg',earthModel);

figure(1),clf
plotpatch2d(boundary,boundary.xc(:,2)), hold on
plotpatch2d(rcv,rcv.xc(:,2))
axis tight equal

%% compute displacement & traction kernels for boundary
[Kdd,Kdn,Knd,Knn] = geometry.computeFullTractionKernels(boundary,boundary);
[Gdx,Gdz,Gnx,Gnz] = computeDisplacementKernels(boundary,boundary.xc);
