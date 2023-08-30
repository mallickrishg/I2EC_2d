%                 Investigate the
%         viscous reponse of oceanic lithosphere 
%          using imposed earthquake cycle method    
%            S. Sathiakumar & Rishav Mallick, 
%               EOS, Caltech, August 2023

clear  
close all
addpath functions/
addpath meshing/
import geometry.*

% Elastic parameters (homogenous medium)
nu=0.25;% Poisson's ratio
mu=30e3;% in MPa

%% geometry
earthModel = geometry.LDhs(mu,nu);

% megathrust fault
Nfault = 50;
Fault_x = 200e3;
Fault_z = 50e3;
x0 = 0;
z0 = 0;
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
