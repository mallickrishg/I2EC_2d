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

%% load fault and shear zone meshes
earthModel = geometry.LDhs(mu,nu);

% megathrust fault
Nfault = 50;
Fault_x = 200e3;
Fault_z = 50e3;
x0 = 0;
z0 = 0.01e3;
dip = atand(Fault_z/Fault_x); % 
Fault_width = sqrt(Fault_x^2 + Fault_z^2);

rcv = create_megathrust(earthModel,x0,z0,dip,Fault_width,Nfault);
shz = load_viscous_wedges(earthModel);

figure(1),clf
plotpatch2d(rcv,rcv.xc(:,2)./1e3), hold on
plotshz2d(shz,shz.xc(:,2)./1e3)
axis tight equal
box on
set(gca,'YDir','normal','Fontsize',20,'Linewidth',2)

%% compute stress interaction kernels

evl = compute_all_stresskernels(rcv,shz);

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





