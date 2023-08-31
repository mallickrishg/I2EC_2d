%                 Investigate the
%         viscous reponse of oceanic lithosphere 
%          using imposed earthquake cycle method    
%            S. Sathiakumar & Rishav Mallick, 
%               EOS, Caltech, August 2023

clear  
addpath functions/
import geometry.*

% Elastic parameters (homogenous medium)
nu=0.25;% Poisson's ratio
mu=30e3;% in MPa

% create megathrust (1) or load existing file (0)
create_fault = 0;

% Periodic earthquake recurrence time
Trecur = 100*3.15e7;% in seconds
Vpl = 1e-9;% m/s
%% load fault and shear zone meshes
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

% provide shear zone mesh as 2 .dat files of the form
% meshname_vertices.dat (contains x,z coordinates of vertices)
% meshname_triangulation.dat (contains 3 columns of vertex linkage)
shz = load_viscous_wedges('inputs/shearzone',earthModel);

figure(1),clf
plotpatch2d(rcv,rcv.xc(:,2)./1e3), hold on
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
evl = compute_all_stresskernels(rcv,shz);

%% assign rheological properties 
% (assuming spatially constant values)
rcv.Asigma = 0.5.*ones(rcv.N,1);% (a-b)sigma
shz.alpha = 1/(1e18*1e-6).*ones(shz.N,1); % alpha = 1/viscosity where viscosity is in MPa-s
shz.n = ones(shz.N,1);

% define locked zone on megathrust
locked = abs(rcv.xc(:,2)) > 15e3 & abs(rcv.xc(:,2))< 40e3;
rcv.pinnedPosition = false(shz.N,1);
rcv.pinnedPosition(locked) = true;

% define long-term slip/strain rates
rcv.Vpl = Vpl;% m/s
shz.e22pl = 1e-15;% 1/s
shz.e23pl = 1e-14;% 1/s

%% calculate coseismic stress change - imposed periodically
slip_coseismic = zeros(rcv.N,1);
slip_coseismic(rcv.pinnedPosition) = Trecur*Vpl;% in meters

delta_stress_rcv = evl.KK*slip_coseismic;
delta_stress22_shz = evl.KL(:,:,1)*slip_coseismic;
delta_stress23_shz = evl.KL(:,:,2)*slip_coseismic;

% plot stress change
figure(2),clf
subplot(3,1,1)
plot(rcv.xc(:,1)./1e3,delta_stress_rcv,'LineWidth',2)
grid on
xlim([-100 350])
xlabel('x (km)'), ylabel('\Delta\tau (MPa)')
subplot(3,1,2)
plotpatch2d(rcv)
plotshz2d(shz,delta_stress22_shz)
axis tight equal
cb=colorbar;cb.Label.String = '\sigma_{xx}^{dev} (MPa)';
clim([-1 1]*0.5)
subplot(3,1,3)
plotpatch2d(rcv)
plotshz2d(shz,delta_stress23_shz)
axis tight equal
cb=colorbar;cb.Label.String = '\sigma_{xz} (MPa)';
clim([-1 1]*0.5)
colormap("bluewhitered")

%% testing dummy
figure(11),clf
toplot = zeros(shz.N,1);
x0 = 0e3;
z0 = -30e3;
r = 15e3;
shzindex = sqrt((shz.xc(:,1)-x0).^2 + (shz.xc(:,2)-z0).^2)<=r;
toplot(shzindex) = 1e-4;

subplot(2,1,1)
plot(rcv.xc(:,1)./1e3,evl.LK(:,:,1)*toplot,'-','Linewidth',2)
axis tight
xlim([-100 350])

subplot(2,1,2)
plotpatch2d(rcv,evl.LK(:,:,1)*toplot)
plotshz2d(shz,evl.LL(:,:,1,1)*toplot)
axis tight equal
box on
% colorbar
clim([-1 1]*max(abs(get(gca,'CLim'))))
% colormap("jet")
colormap("bluewhitered")





