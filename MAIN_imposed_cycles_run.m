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
Nfault = 100;
Fault_width = 200e3;
Fault_depth = 50e3;
x0 = 0;
z0 = 0;
dip = asind(Fault_depth/Fault_width); % 

rcv = create_megathrust(earthModel,x0,z0,dip,Fault_width,Nfault);
shz = load_viscous_wedges('viscous_shz');




