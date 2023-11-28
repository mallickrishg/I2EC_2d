%% Plot the long term velocity field

tic
clear  
addpath functions/
addpath inputs/
import geometry.*


% Elastic parameters (homogenous medium)
nu=0.25;% Poisson's ratio
mu=30e3;% in MPa

% create megathrust (1) or load existing file (0)
create_fault = 0;

Vpl = 1e-9;% m/s
Tplate=20e3;

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
boundary = geometry.receiver('inputs/boundary2d.seg',earthModel);
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
plotshz2d(shz)
axis tight equal
box on
set(gca,'YDir','normal','Fontsize',20,'Linewidth',2)

%% Observation points
x=linspace(min( boundary.x(:,1)),max( boundary.x(:,1)),50);
y=linspace(max( boundary.x(:,2)),min( boundary.x(:,2)),50);
[X,Y]=meshgrid(x,y);
X=X(:); Y=Y(:); dip=rcv.dip(1);
ocean_translate_params = [0,20e3];
continental_translate_params=[-140e3,35e3];

[Vx,Vz]= computeLongTermVelocity(X,Y,Vpl,Tplate,dip,shz,boundary,ocean_translate_params,continental_translate_params);


hold on;
axis equal
plot(X/1e3,Y/1e3,'c.')
quiver(X/1e3,Y/1e3,Vx,Vz,'k')
