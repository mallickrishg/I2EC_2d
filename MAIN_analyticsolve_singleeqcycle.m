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
Trecur = 500*3.15e7;% in seconds
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
% evl_orig = computeAllStressKernelsBem(rcv,shz,boundary,'kernelmodify',0);
load('kernels/evl_orig.mat','evl_orig');

% compute displacement kernels
Nobs = 401;
obs = ([1;0]*(linspace(-100,500,Nobs)))'*1e3;
% devl = computeAllDisplacementKernelsBem(obs,rcv,shz,boundary,1);
load('kernels/devl.mat','devl');

%% assign rheological properties 

%%%%%%% approximate fault by a viscous shear zone %%%%%%%
% (a-b)sigma in terms of viscosity
% eta' = viscosity/L_fault
rcv.Asigma = 1e-6.*(1e18/(sum(~rcv.pinnedPosition.*rcv.W))).*ones(rcv.N,1);

% impose viscosity structure
shz.n = 1.*ones(shz.N,1);
%%%%%%% oceanic mantle viscosity structure %%%%%%%
r = abs(tand(rcv.dip(1)).*shz.xc(:,1) + shz.xc(:,2) + 20e3)./sqrt(tand(rcv.dip(1))^2 + 1);
r = r./max(r);% normalize to 0->1
viscostructure = 10.^(18 + r.*0);
shz.alpha = 1./(viscostructure.*1e-6);

oceanic_mantle = (shz.xc(:,1) < -shz.xc(:,2)/tand(rcv.dip(1)));
% %%%%%%% continental mantle viscosity structure %%%%%%%
r = sqrt((shz.xc(~oceanic_mantle,1)-200e3).^2);
r = r./max(r);% normalize to 0->1
viscostructure = 10.^(18 + (((r-0.5).*2).^2).*0);
shz.alpha(~oceanic_mantle) = 1./(viscostructure.*1e-6);

% define locked zone on megathrust
locked = abs(rcv.xc(:,2)) > 0e3 & abs(rcv.xc(:,2))< 30e3;
rcv.pinnedPosition = locked;

% define long-term slip/strain rates
rcv.Vpl = Vpl.*ones(rcv.N,1);% m/s

% Long-term strain rate calculation
[e22_dev_lt, e23_lt] = getStrainratesLongterm(shz,rcv.dip(1)*pi/180,[0,20e3],[-140e3,35e3],1);
shz.e22pl = e22_dev_lt.*Vpl;% 1/s
shz.e23pl = -e23_lt.*Vpl;% 1/s

figure(1),clf
plotshz2d(shz,1e6./(shz.alpha))
axis tight equal, box on
cb=colorbar; cb.Label.String = '\eta Pa-s';
set(gca,'ColorScale','log','FontSize',15,'LineWidth',1.5,'TickDir','both')

%% calculate coseismic stress change - imposed periodically
slip_coseismic = zeros(rcv.N,1) + Trecur*Vpl.*(locked);

% initialise stress change data structure
stress_change = [];

stress_change.dtau = zeros(rcv.N,1);
stress_change.dsigma22 = zeros(shz.N,1);
stress_change.dsigma23 = zeros(shz.N,1);

% stress change for the coseismic event 
dtau = evl_orig.KK*slip_coseismic;
stress_change.dtau = dtau;
stress_change.dsigma22 = evl_orig.KL(:,:,1)*slip_coseismic;
stress_change.dsigma23 = evl_orig.KL(:,:,2)*slip_coseismic;

%% construct stress interaction kernel (use unmodified stress kernel for interseismic calculation)
stresskernel = assembleStressKernel(evl_orig,locked);
viscosityvector = [rcv.Asigma(~locked);1./shz.alpha;1./shz.alpha];
% construct Operator matrix for dy/dt = [A].y 
% where [A] = [K/η] - linear operator that contains all rheological parameters
rheoparam = diag(1./viscosityvector)*stresskernel;

% combine all stress change into a single vector that captures coseismic
% strain rate change
deltastrainrate = [stress_change.dtau(~locked);...
                   stress_change.dsigma22;...
                   stress_change.dsigma23]./viscosityvector;

% compute initial condition (and late-interseismic condition)
fault_locking_velocity = zeros(rcv.N,1);
fault_locking_velocity(locked) = Vpl;
% strain rate loading from locking + long-term
edot_locking = -[evl_orig.KK(~locked,:)*fault_locking_velocity;...
          evl_orig.KL(:,:,1)*fault_locking_velocity;...
          evl_orig.KL(:,:,2)*fault_locking_velocity]./viscosityvector;
edot_lt = -rheoparam*[rcv.Vpl(~locked);shz.e22pl;shz.e23pl];

% initial condition
sol_initial = computePeriodicInitialCondition(rheoparam, edot_locking + edot_lt,...
                                         deltastrainrate, Trecur);
% late-interseismic condition
sol_interseismic = sol_initial - deltastrainrate;

%% provide time-steps to evaluate the solution
tvec = [40,365,10*365,365*50].*86400;

% use reconditioned matrix if singularity issues
[sol, sol_integrated] = analyticSolveIVP(rheoparam, edot_locking + edot_lt,...
    sol_initial, tvec);

[V,e22dot,e23dot] = extractComponentsSolutionVector(sol,locked,rcv.N,shz.N);
[slip,e22,e23] = extractComponentsSolutionVector(sol_integrated,locked,rcv.N,shz.N);

%% plot snapshots of internal strain rate
figure(2),clf
for i = 1:length(tvec)
    tval = tvec(i);
    % subplot(length(tvec)-1,1,i)
    subplot(2,2,i)
    plotshz2d(shz,sqrt(e22dot(:,i).^2 + e23dot(:,i).^2))
    axis tight equal, box on
    clim(10.^[-15 -12])
    cb=colorbar; cb.Label.String = 'strain rate 1/s';
    colormap(turbo(30))
    title(['\Deltat = ' num2str(round(tval/86400/365,2)) ' years'])
    set(gca,'ColorScale','log','FontSize',12,'LineWidth',1.5,'TickDir','in')
end

% plot internal strain rates
figure(3),clf
% long-term
subplot(3,1,1)
plotshz2d(shz,sqrt(shz.e22pl.^2 + shz.e23pl.^2))
axis tight equal, box on
clim(10.^[-15 -13])
cb=colorbar; cb.Label.String = 'strain rate (1/s)';
xlabel('x (km)'),ylabel('depth (km)')
title('$\dot{\epsilon}_v^{\infty}$ ($t \to \infty$)','Interpreter','latex')
set(gca,'ColorScale','log','FontSize',15,'LineWidth',1,'TickDir','in')

% initial-condition
e22_0 = sol_initial(length(find(~locked))+1:length(find(~locked))+shz.N);
e23_0 = sol_initial(length(find(~locked))+shz.N+1:end);
subplot(3,1,2)
plotshz2d(shz,sqrt(e22_0.^2 + e23_0.^2))
axis tight equal, box on
clim(10.^[-14 -12])
cb=colorbar; cb.Label.String = 'strain rate (1/s)';
xlabel('x (km)'),ylabel('depth (km)')
title('$\dot{\epsilon}_v$ ($t = 0$)','Interpreter','latex')
set(gca,'ColorScale','log','FontSize',15,'LineWidth',1,'TickDir','in')

% late-interseismic
e22_T = sol_interseismic(length(find(~locked))+1:length(find(~locked))+shz.N);
e23_T = sol_interseismic(length(find(~locked))+shz.N+1:end);
subplot(3,1,3)
plotshz2d(shz,sqrt(e22_T.^2 + e23_T.^2))
axis tight equal, box on
clim(10.^[-15 -13])
cb=colorbar; cb.Label.String = 'strain rate (1/s)';
colormap("hot")
xlabel('x (km)'),ylabel('depth (km)')
title('$\dot{\epsilon}_v$ ($t = T_{eq}$)','Interpreter','latex')
set(gca,'ColorScale','log','FontSize',15,'LineWidth',1,'TickDir','in')

% plot late interseismic velocity on the fault
figure(4),clf
plot(rcv.xc(~locked,1)./1e3,sol_interseismic(1:length(find(~locked)))./Vpl,'-','LineWidth',2)
axis tight, 
xlabel('x (km)'), ylabel('v/v_{pl}')
xlim([0 max(rcv.xc(:,1)./1e3)])
ylim([0 1])

% return
%% surface interseismic velocities
% calculate displacement time series at select observation points
stationID = (100:50:350);

fullsol_interseismic = [zeros(length(find(locked)),1);sol_interseismic] - [rcv.Vpl;shz.e22pl;shz.e23pl];
gps = [];
gps.vx_int = [devl.KO(:,:,1),devl.LO(:,:,1,1),devl.LO(:,:,1,2)]*fullsol_interseismic - Ghingex_d * (hinge.Vpl.*Vpl); 
gps.vz_int = [devl.KO(:,:,2),devl.LO(:,:,2,1),devl.LO(:,:,2,2)]*fullsol_interseismic - Ghingez_d * (hinge.Vpl.*Vpl); 
% displacements
gps.ux = (devl.KO(:,:,1)*(slip-rcv.Vpl*tvec) + ...
          devl.LO(:,:,1,1)*(e22-shz.e22pl*tvec) + ... 
          devl.LO(:,:,1,2)*(e23-shz.e23pl*tvec) - ...
          1.*Ghingex_d * (hinge.Vpl.*Vpl)*tvec)';
gps.uz = (devl.KO(:,:,2)*(slip-rcv.Vpl*tvec) + ...
          devl.LO(:,:,2,1)*(e22-shz.e22pl*tvec) + ...
          devl.LO(:,:,2,2)*(e23-shz.e23pl*tvec) - ...
          1.*Ghingez_d * (hinge.Vpl.*Vpl)*tvec)';
% interseismic velocity
figure(12),clf
plot(obs(:,1)./1e3,gps.vx_int./Vpl,'LineWidth',2), hold on
plot(obs(:,1)./1e3,gps.vz_int./Vpl,'LineWidth',2)
axis tight,grid on
ylim([-1 1]*1)
legend('horizontal','vertical','Box','off','Location','best')
xlabel('distance from trench (km)'), ylabel('v/v_{pl}')
set(gca,'FontSize',15,'LineWidth', 1,'TickDir','both')

% displacement snapshots
figure(13),clf
subplot(2,1,1)
plot(obs(:,1)./1e3,gps.ux,'LineWidth',2), hold on
plot(obs(stationID,1)./1e3,zeros(length(stationID),1),'ko','MarkerFaceColor','k')
axis tight,grid on
% ylim([-1 1]*1)
xlabel('distance from trench (km)'), ylabel('u_h (m)')
set(gca,'FontSize',15,'LineWidth', 1,'TickDir','both')
subplot(2,1,2)
plot(obs(:,1)./1e3,gps.uz,'LineWidth',2), hold on
plot(obs(stationID,1)./1e3,zeros(length(stationID),1),'ko','MarkerFaceColor','k')
axis tight,grid on
% ylim([-1 1]*1)
xlabel('distance from trench (km)'), ylabel('u_z (m)')
set(gca,'FontSize',15,'LineWidth', 1,'TickDir','both')
% return
%% plot time series at selected stations
cspec = parula(length(stationID));
npts = 40;
t_inter = Trecur - linspace(5,0,npts).*3.15e7;
t_post = logspace(5,8.3,npts);
tvec = [t_inter,t_post];

% solve for slip & strain time series
tic
[~, sol_integrated] = analyticSolveIVP(rheoparam, edot_locking + edot_lt,...
    sol_initial, tvec);
[slip,e22,e23] = extractComponentsSolutionVector(sol_integrated,locked,rcv.N,shz.N);
toc
% displacements
gps.ux = (devl.KO(:,:,1)*(slip-rcv.Vpl*tvec) + ...
          devl.LO(:,:,1,1)*(e22-shz.e22pl*tvec) + ... 
          devl.LO(:,:,1,2)*(e23-shz.e23pl*tvec) - ...
          1.*Ghingex_d * (hinge.Vpl.*Vpl)*tvec)';
gps.uz = (devl.KO(:,:,2)*(slip-rcv.Vpl*tvec) + ...
          devl.LO(:,:,2,1)*(e22-shz.e22pl*tvec) + ...
          devl.LO(:,:,2,2)*(e23-shz.e23pl*tvec) - ...
          1.*Ghingez_d * (hinge.Vpl.*Vpl)*tvec)';


plotshift = 0.5;
figure(14),clf
subplot(2,1,1)
for i = 1:length(stationID)
    toplot = gps.ux(1:npts,stationID(i));
    plot((t_inter-Trecur)./3.15e7,toplot-toplot(end)+(i-1)*plotshift,'.-','LineWidth',2,'Color',cspec(i,:))
    hold on
    toplot = gps.ux(1+npts:end,stationID(i));
    plot(t_post./3.15e7,toplot+(i-1)*plotshift,'.-','LineWidth',2,'Color',cspec(i,:))
    axis tight
end
xlabel('time (yrs)')
ylabel('u_h (m)')
set(gca,'FontSize',15,'LineWidth',1)
plotshift = 1;
subplot(2,1,2)
for i = 1:length(stationID)
    toplot = gps.uz(1:npts,stationID(i));
    plot((t_inter-Trecur)./3.15e7,toplot-toplot(end)+(i-1)*plotshift,'.-','LineWidth',2,'Color',cspec(i,:))
    hold on
    toplot = gps.uz(1+npts:end,stationID(i));
    plot(t_post./3.15e7,toplot+(i-1)*plotshift,'.-','LineWidth',2,'Color',cspec(i,:))
    axis tight
end
xlabel('time (yrs)')
ylabel('u_z (m)')
set(gca,'FontSize',15,'LineWidth',1)

