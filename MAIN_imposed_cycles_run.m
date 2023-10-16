%                 Investigate the
%         viscous reponse of oceanic lithosphere 
%          using imposed earthquake cycle method    
%            S. Sathiakumar & Rishav Mallick, 
%               EOS, Caltech, August 2023

clear  
addpath functions/
import geometry.*
rng(42)

% Elastic parameters (homogenous medium)
nu=0.25;% Poisson's ratio
mu=30e3;% in MPa

% create megathrust (1) or load existing file (0)
create_fault = 0;

% Periodic earthquake recurrence time
Trecur = 100*3.15e7;% in seconds
Vpl = 1e-9;% m/s

% max stress change on fault (MPa)
tau_max = 5;
%% load fault, boundary and shear zone meshes
earthModel = geometry.LDhs(mu,nu);

% megathrust fault
if create_fault == 1
    Nfault = 50;
    Fault_x = 200e3;
    Fault_z = 50e3;
    x0 = 0;
    z0 = 0.01e3;
    dip = atand(Fault_z/Fault_x); %
    Fault_width = sqrt(Fault_x^2 + Fault_z^2);
    % create megathrust mesh
    rcv = create_megathrust(earthModel,x0,z0,dip,Fault_width,Nfault);
else
    % directly load existing file
    rcv = geometry.receiver('inputs/megathrust2d.seg',earthModel);
end

% boundary mesh
boundary = geometry.receiver('inputs/boundarytight2d.seg',earthModel);
boundary.Vx = boundary.Vx.*Vpl;
boundary.Vz = boundary.Vz.*Vpl;

% provide shear zone mesh as 2 .dat files of the form
% meshname_vertices.dat (contains x,z coordinates of vertices)
% meshname_triangulation.dat (contains 3 columns of vertex linkage)
% This mesh can be created using CREATE_shearzone_mesh.m provided in the
% folder 'meshing'
shz = geometry.shearZoneReceiver('inputs/shearzone',earthModel);

figure(1),clf
plotpatch2d(rcv,rcv.xc(:,2)./1e3), hold on
plotpatch2d(boundary,boundary.xc(:,2)./1e3)
plotshz2d(shz,shz.xc(:,2)./1e3)
axis tight equal
box on
set(gca,'YDir','normal','Fontsize',20,'Linewidth',2)

%% compute stress interaction kernels
% evl contains the following as N-d matrices
% KK - fault-fault interactions [rcv.N x rcv.N]
% KL - fault-shz interactions [shz.N x rcv.N x 2]
% LK - shz-fault interactions [rcv.N x shz.N x 2]
% LL - shz-shz interactions [shz.N x shz.N x 2 x 2]
evl = computeAllStressKernelsBem(rcv,shz,boundary);

%% assign rheological properties 
% (assuming spatially constant values)
rcv.Asigma = 0.5.*ones(rcv.N,1);% (a-b)sigma
shz.alpha = 1/(2e20*1e-6).*ones(shz.N,1); % alpha = 1/viscosity where viscosity is in MPa-s
shz.n = ones(shz.N,1)+0.1;

% define locked zone on megathrust
locked = abs(rcv.xc(:,2)) > 10e3 & abs(rcv.xc(:,2))< 40e3;
rcv.pinnedPosition = false(rcv.N,1);
rcv.pinnedPosition(locked) = true;

% define long-term slip/strain rates
rcv.Vpl = Vpl.*ones(rcv.N,1);% m/s

% Long-term strain rate calculation
[e22_dev, e23] = getStrainratesLongterm(shz,rcv.dip(1)*pi/180,[0,20e3],[-140e3,35e3]);
shz.e22pl = e22_dev.*Vpl;
shz.e23pl = e23.*Vpl;
% shz.e22pl = 1e-15.*ones(shz.N,1);% 1/s
% shz.e23pl = 1e-14.*ones(shz.N,1);% 1/s

%% calculate coseismic stress change - imposed periodically
Nevents = 1;
slip_coseismic = zeros(rcv.N,Nevents);
slip_multipliers = drchrnd(ones(1,Nevents),1);% this is just to create random numbers that sum to 1

for i = 1:Nevents
    slip_coseismic(rcv.pinnedPosition,i) = Trecur*Vpl.*slip_multipliers(i);% in meters
end

% initialise stress change data structure
stress_change = [];
stress_change.Nevents = Nevents;
stress_change.Timing = 4*3.15e7;%[4,10,50]*3.15e7;% provide earthquake timing (in seconds) as a vector

assert(length(stress_change.Timing) == Nevents)

stress_change.dtau = zeros(rcv.N,stress_change.Nevents);
stress_change.dsigma22 = zeros(shz.N,stress_change.Nevents);
stress_change.dsigma23 = zeros(shz.N,stress_change.Nevents);

% stress change for each event stored as a matrix
for i = 1:stress_change.Nevents
    dtau = evl.KK*slip_coseismic(:,i);
    dtau(dtau > tau_max) = tau_max;
    stress_change.dtau(:,i) = dtau;
    stress_change.dtau(locked,i) = 0;% force stress change in coseismic region to 0

    stress_change.dsigma22(:,i) = evl.KL(:,:,1)*slip_coseismic(:,i);
    stress_change.dsigma23(:,i) = evl.KL(:,:,2)*slip_coseismic(:,i);
end

% plot stress change
figure(2),clf
subplot(3,1,1)
plot(rcv.xc(:,1)./1e3,stress_change.dtau(:,1),'LineWidth',2)
grid on
xlim([-100 350])
xlabel('x (km)'), ylabel('\Delta\tau (MPa)')
subplot(3,1,2)
plotpatch2d(rcv)
plotshz2d(shz,stress_change.dsigma22(:,1))
axis tight equal
cb=colorbar;cb.Label.String = '\sigma_{xx}^{dev} (MPa)';
clim([-1 1]*0.5)
subplot(3,1,3)
plotpatch2d(rcv)
plotshz2d(shz,stress_change.dsigma23(:,1))
axis tight equal
cb=colorbar;cb.Label.String = '\sigma_{xz} (MPa)';
clim([-1 1]*0.5)
colormap("bluewhitered")

% return
%% use rcv, evl, shz, stress_change to run earthquake cycles
Ncycles = 5;% specify number of cycles (for spin up)
tic
disp('running imposed earthquake sequence simulations')
[t,V,e22dot,e23dot] = runImposedEarthquakeCycles(rcv,shz,evl,stress_change,Ncycles,Trecur);
toc
%% plot results
edot_pl = sqrt(shz.e22pl.^2 + shz.e23pl.^2);
edot_pl = mean(edot_pl).*ones(shz.N,1);
edot = sqrt(e22dot.^2 + e23dot.^2);

figure(10),clf
pcolor(t./Trecur,rcv.xc(:,1)./1e3,V'./Vpl), shading interp
xlabel('t/T_{eq}')
ylabel('x (km)')
colorbar
clim(10.^[-1,2])
colormap("turbo")
set(gca,'ColorScale','log','YDir','reverse','FontSize',15,'TickDir','out','LineWidth',1.5)

figure(11),clf
shzindex = find(sqrt((shz.xc(:,1)-50e3).^2 + (shz.xc(:,2)+40e3).^2) < 10e3);
plot(t./Trecur,V(:,43)./Vpl,'.-'), hold on
for i = 1:length(shzindex)
    plot(t./Trecur,edot(:,shzindex(i))./edot_pl(shzindex(i)),'r.-')
end
axis tight
legend('slip rate','strain rate')
xlabel('t/T_{eq}')
set(gca,'YScale','log','FontSize',15,'TickDir','out','LineWidth',1.5)
ylabel('$\frac{v}{v_{pl}}$ , $\frac{\dot{\epsilon}}{\dot{\epsilon}_{pl}}$','Interpreter','latex','FontSize',25)

figure(100),clf
toplot = zeros(shz.N,1);
toplot(shzindex) = 1;
plotshz2d(shz,toplot)
axis tight equal
%% create snapshots of normalized slip rate & strain rates
% t_plots = [0,4.5,5.01,6,10,19.5].*3.15e7;
% plotindex = [5,9,11,49,51,99].*3.15e7;
plotindex = [7,10,15,20,50,80].*3.15e7;
figure(12),clf
for i = 1:length(plotindex)
    tindex = find(abs(t-plotindex(i))==min(abs(t-plotindex(i))),1);
    subplot(3,2,i)
    plotshz2d(shz,edot(tindex,:)'./edot_pl), hold on
    plotpatch2d(rcv,V(tindex,:)'./rcv.Vpl)
    box on
    cb=colorbar;cb.Label.String = '\gamma/\gamma_0';
    colormap("turbo")
    clim(10.^[-1.5,1.5])
    axis tight equal
    xlabel('x (km)'), ylabel('z (km)')
    title(['t = ' num2str(t(tindex)./3.15e7,'%.1f') ' yrs'])
    set(gca,'Fontsize',15,'ColorScale','log','Linewidth',1.5,'TickDir','out')
end

% return
%% calculate velocity time series at select observation points
Nobs = 1000;
obs = ([1;0]*(linspace(-100,350,Nobs)))'*1e3;

% compute displacement kernels
% devl = compute_all_dispkernels(obs,rcv,shz,boundary,Vpl);

gps = computeSiteVelocitiesBem(obs,rcv,shz,boundary,V,e22dot,e23dot);

%% plotting surface displacements 
figure(3);clf
p = [];
lgd = {};
% plotindex = [0,4,5.01,6,10,19.5].*3.15e7;
cspec = cool(length(plotindex));

subplot(2,1,1);hold on;
toplot=gps.vx;
toplot_pl=Vpl;
plot(obs(:,1)./1e3,(toplot(end,:))./toplot_pl,'k-','LineWidth',3), hold on
for i = 1:length(plotindex)
    tindex = find(abs(t-plotindex(i))==min(abs(t-plotindex(i))),1);
    plot(obs(:,1)./1e3,(toplot(tindex,:))./toplot_pl,'-','LineWidth',2,'Color',cspec(i,:));
end
axis tight
grid on;box on
xlabel('distance from trench (km)'), ylabel('v_x/v_{pl}')
title("Horizontal component")

subplot(2,1,2); hold on;
toplot=gps.vz;
toplot_pl=Vpl;
plot(obs(:,1)./1e3,(toplot(end,:))./toplot_pl,'k-','LineWidth',3), hold on
for i = 1:length(plotindex)
    tindex = find(abs(t-plotindex(i))==min(abs(t-plotindex(i))),1);
    p(i) = plot(obs(:,1)./1e3,(toplot(tindex,:))./toplot_pl,'-','LineWidth',2,'Color',cspec(i,:));
    lgd{i} = [num2str(round(plotindex(i)./3.15e7)) ' yrs'];
end
legend(p,lgd); 
axis tight
grid on;box on
xlabel('distance from trench (km)'), ylabel('v_z/v_{pl}')
title("Vertical component")

set(findobj(gcf,'type','axes'),'FontSize',15,'LineWidth', 1);

%% plot velocity cross-sections as snapshots
x = linspace(-100,349,40).*1e3;
z = linspace(-79,0,10).*1e3;
[X,Z] = meshgrid(x,z);
obs = [X(:),Z(:)];
gps = computeSiteVelocitiesBem(obs,rcv,shz,boundary,V,e22dot,e23dot);

figure(13),clf
for i = 1:length(plotindex)
    tindex = find(abs(t-plotindex(i))==min(abs(t-plotindex(i))),1);
    vtot = reshape(sqrt(gps.vx(tindex,:).^2 + gps.vz(tindex,:).^2),length(z),length(x));
    
    subplot(3,2,i)
    pcolor(x./1e3,z./1e3,vtot./Vpl), shading interp
    alpha(0.5)
    hold on
    plotshz2d(shz)
    plotpatch2d(rcv)
    quiver(obs(:,1)./1e3,obs(:,2)./1e3,gps.vx(tindex,:)'./Vpl,gps.vz(tindex,:)'./Vpl,'k','LineWidth',1)
    axis tight equal
    cb=colorbar;cb.Label.String = 'v/v_{pl}';
    colormap("turbo")
    clim([0 1.2])
    colormap(parula(12))
    xlabel('x (km)'), ylabel('z (km)')
    title(['t = ' num2str(t(tindex)./3.15e7,'%.1f') ' yrs'])
    set(gca,'Fontsize',15,'ColorScale','lin','Linewidth',1.5,'TickDir','out')
end

figure(14),clf
for i = 1:length(plotindex)
    tindex = find(abs(t-plotindex(i))==min(abs(t-plotindex(i))),1);
    vtot = reshape(gps.vz(tindex,:),length(z),length(x));
    
    subplot(3,2,i)
    pcolor(x./1e3,z./1e3,vtot./Vpl), shading interp
    hold on
    plotshz2d(shz)
    plotpatch2d(rcv)
    quiver(obs(:,1)./1e3,obs(:,2)./1e3,gps.vx(tindex,:)'./Vpl,gps.vz(tindex,:)'./Vpl,'k','LineWidth',1)
    axis tight equal
    cb=colorbar;cb.Label.String = 'v_z/v_{pl}';
    clim([-1 1])
    colormap(bluewhitered(20))
    xlabel('x (km)'), ylabel('z (km)')
    title(['t = ' num2str(t(tindex)./3.15e7,'%.1f') ' yrs'])
    set(gca,'Fontsize',15,'ColorScale','lin','Linewidth',1.5,'TickDir','out')
end
%% steady state motion

gps_ss = computeSiteVelocitiesBem(obs,rcv,shz,boundary,rcv.Vpl',shz.e22pl',shz.e23pl');
%%
toplot = reshape(gps_ss.vx(1,:),length(z),length(x))./Vpl;

figure(4),clf
plot(x./1e3,toplot(10,:),'-','Linewidth',2)
axis tight, grid on