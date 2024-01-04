% Main script for 
% Earthquake cycle calculation using 
% semi-analytical solutions for a set of coupled ODEs
% 
% AUTHOR:
% Rishav Mallick, JPL 2024

clear  
addpath functions/
import('geometry.*')

% Elastic parameters (homogenous medium)
nu=0.25;% Poisson's ratio
mu=30e3;% in MPa

% Periodic earthquake recurrence time
Trecur = 100*3.15e7;% in seconds
Vpl = 1e-9;% m/s

% max stress change on fault (MPa)
tau_max = 3;
%% load fault, boundary and shear zone meshes
earthModel = geometry.LDhs(mu,nu);

rcv = geometry.receiver('inputs/megathrust2d.seg',earthModel);

% boundary mesh
boundary = geometry.receiver('inputs/boundary2d.seg',earthModel);
boundary.Vx = boundary.Vx.*Vpl;
boundary.Vz = boundary.Vz.*Vpl;

% provide shear zone mesh as 2 .dat files of the form
% meshname_vertices.dat (contains x,z coordinates of vertices)
% meshname_triangulation.dat (contains 3 columns of vertex linkage)
% This mesh can be created using CREATE_shearzone_mesh.m provided in the
% folder 'meshing'

shz = geometry.shearZoneReceiver('inputs/shearzone',earthModel);

%% compute stress interaction kernels
% evl contains the following as N-d matrices
% KK - fault-fault interactions [rcv.N x rcv.N]
% KL - fault-shz interactions [shz.N x rcv.N x 2]
% LK - shz-fault interactions [rcv.N x shz.N x 2]
% LL - shz-shz interactions [shz.N x shz.N x 2 x 2]

% use original unmodified kernels for this solve
evl_orig = computeAllStressKernelsBem(rcv,shz,boundary,'kernelmodify',0);
% load('kernels/evl_orig.mat','evl_orig');

% compute displacement kernels
Nobs = 401;
obs = ([1;0]*(linspace(-100,500,Nobs)))'*1e3;
devl = computeAllDisplacementKernelsBem(obs,rcv,shz,boundary,1);
% load('kernels/devl.mat','devl');

%% assign rheological properties 

%%%%%%% approximate fault by a viscous shear zone %%%%%%%
% (a-b)sigma in terms of viscosity
% eta' = viscosity/L_fault
rcv.Asigma = 1e-6.*(1e19/(sum(~rcv.pinnedPosition.*rcv.W))).*ones(rcv.N,1);

%%%%%%% oceanic mantle viscosity structure %%%%%%%
r = abs(tand(rcv.dip(1)).*shz.xc(:,1) + shz.xc(:,2) + 20e3)./sqrt(tand(rcv.dip(1))^2 + 1);
r = r./max(r);% normalize to 0->1
viscostructure = 10.^(19 + r.*3);
shz.n = 1.*ones(shz.N,1);
shz.alpha = 1./(viscostructure.*1e-6);
oceanic_mantle = (shz.xc(:,1) < -shz.xc(:,2)/tand(rcv.dip(1)));

%%%%%%% continental mantle viscosity structure %%%%%%%
r = sqrt((shz.xc(~oceanic_mantle,1)-200e3).^2);
r = r./max(r);% normalize to 0->1
viscostructure = 10.^(18 + r.*3);
shz.alpha(~oceanic_mantle) = 1./(viscostructure.*1e-6);

% define locked zone on megathrust
locked = abs(rcv.xc(:,2)) > 0e3 & abs(rcv.xc(:,2))< 30e3;
rcv.pinnedPosition = locked;

% define long-term slip/strain rates
rcv.Vpl = Vpl.*ones(rcv.N,1);% m/s

% compute steady state interseismic slip rate & strain rates
[v_ss,e22_ss,e23_ss] = computeInterseismicSteadystate(rcv,shz,Vpl,10);

% Long-term strain rate calculation
[e22_dev_lt, e23_lt] = getStrainratesLongterm(shz,rcv.dip(1)*pi/180,[0,20e3],[-140e3,35e3]);
shz.e22pl = e22_dev_lt.*Vpl;% 1/s
shz.e23pl = -e23_lt.*Vpl;% 1/s

figure(1),clf
plotshz2d(shz,1e6./(shz.alpha))
axis tight equal, box on
cb=colorbar; cb.Label.String = '\eta Pa-s';
set(gca,'ColorScale','log','FontSize',15,'LineWidth',1.5,'TickDir','both')

%% calculate coseismic stress change - imposed periodically
Nevents = 1;
slip_coseismic = zeros(rcv.N,Nevents) + Trecur*Vpl.*(locked);

% initialise stress change data structure
stress_change = [];
stress_change.Nevents = Nevents;
stress_change.Timing = 1;
assert(length(stress_change.Timing) == Nevents)

stress_change.dtau = zeros(rcv.N,stress_change.Nevents);
stress_change.dsigma22 = zeros(shz.N,stress_change.Nevents);
stress_change.dsigma23 = zeros(shz.N,stress_change.Nevents);

% stress change for each event stored as a matrix
for i = 1:stress_change.Nevents
    dtau = evl_orig.KK*slip_coseismic(:,i);
    dtau(dtau > tau_max) = tau_max;
    stress_change.dtau(:,i) = dtau;
    stress_change.dtau(locked,i) = 0;% force stress change in coseismic region to 0

    stress_change.dsigma22(:,i) = evl_orig.KL(:,:,1)*slip_coseismic(:,i);
    stress_change.dsigma23(:,i) = evl_orig.KL(:,:,2)*slip_coseismic(:,i);
end

%% construct stress interaction kernel (use unmodified stress kernel for interseismic calculation)
stresskernel_orig = [evl_orig.KK(~locked,~locked),  evl_orig.LK(~locked,:,1),   evl_orig.LK(~locked,:,2);...
                     evl_orig.KL(:,~locked,1),      evl_orig.LL(:,:,1,1),       evl_orig.LL(:,:,1,2);...
                     evl_orig.KL(:,~locked,2),      evl_orig.LL(:,:,2,1),       evl_orig.LL(:,:,2,2)];

viscosityvector = [rcv.Asigma(~locked);1./shz.alpha;1./shz.alpha];
Nvec = length(viscosityvector);

rheoparam = stresskernel_orig./repmat(viscosityvector,1,Nvec);

% eigenvector decomposition
[Evector,Evals] = eig(rheoparam);
% store eigenvalues in a vector
lambda = diag(Evals);

% combine all stress change into a single vector that captures coseismic
% strain rate change
deltastrainrate = [stress_change.dtau(~locked);stress_change.dsigma22;stress_change.dsigma23]./viscosityvector;

% long-term steady state strain rates
longterm_ss = [v_ss(~locked);e22_ss;e23_ss];

%% compute late interseismic strain rates and initial conditions
% compute late-interseismic strain rate
sol_interseismic = real((eye(Nvec) - Evector*diag(exp(lambda.*Trecur))/Evector)\...
    (Evector*diag(exp(lambda.*Trecur))/Evector*(deltastrainrate-longterm_ss) + longterm_ss));

sol_initial = sol_interseismic + deltastrainrate;

%% evaluate slip/strain at various times

tvec = [7,30,180,365,2*365,5*365,10*365].*86400;% seconds

slip = zeros(rcv.N,length(tvec));
e22 = zeros(shz.N,length(tvec));
e23 = zeros(shz.N,length(tvec));

for i = 1:length(tvec)
    tval = tvec(i);
    sol = real(Evector*diag((exp(lambda.*tval) - ones(Nvec,1))./lambda)/Evector*(sol_initial-sol_interseismic)) + ...
               sol_interseismic.*tval;
    
    % extract solution to slip & strain components
    slip(~locked,i) = sol(1:length(find(~locked)));
    e22(:,i) = sol(length(find(~locked))+1:length(find(~locked))+shz.N);
    e23(:,i) = sol(length(find(~locked))+shz.N+1:end);
end

%% surface displacements
% calculate velocity time series at select observation points (interseismic removed)

hinge = geometry.receiver('inputs/hinge2d.seg',earthModel);
[Gx_d,Gz_d] = computeFaultDisplacementKernelsBem(hinge,obs,boundary,1);

gps = [];
gps.obs = obs;
gps.ux = (devl.KO(:,:,1)*(slip-rcv.Vpl*tvec) + ...
          devl.LO(:,:,1,1)*(e22-shz.e22pl*tvec) + ... 
          devl.LO(:,:,1,2)*(e23-shz.e23pl*tvec) - ...
          1.*Gx_d * (hinge.Vpl.*Vpl)*tvec)';
gps.uz = (devl.KO(:,:,2)*(slip-rcv.Vpl*tvec) + ...
          devl.LO(:,:,2,1)*(e22-shz.e22pl*tvec) + ...
          devl.LO(:,:,2,2)*(e23-shz.e23pl*tvec) - ...
          1.*Gz_d * (hinge.Vpl.*Vpl)*tvec)';

figure(3),clf
subplot(2,1,1)
plot(obs(:,1)./1e3,gps.ux)
axis tight

subplot(2,1,2)
plot(obs(:,1)./1e3,gps.uz)
axis tight


