%                 Investigate the
%         viscous reponse of oceanic lithosphere 
%          using imposed earthquake cycle method    
%            S. Sathiakumar, EOS, August 2023
%%
clear all; close all;
addpath '/Users/sharadha/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/projects/postseismic/matlab/unicycle/matlab';
addpath '/Users/sharadha/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/projects/postseismic/matlab/greens';
import unicycle.*
%%
% Poisson's ratio
nu=0.25;
rig=30e3;
%% geometry

% megathrust fault 
% input fault parameters 

flt.N=100; % no of fault patches 
flt.W = 250e3*ones(flt.N,1); % width of fault 
flt.dip=10*ones(flt.N,1); %dip of fault 
A=[  0;    0]; % starting point 
B=[cosd(flt.dip(1))*flt.W(1); sind(flt.dip(1))*flt.W(1);]; % ending point 

% fault patches 
flt.x  =repmat((0:(flt.N-1))'/flt.N,1,2).*repmat((B-A)',flt.N,1);
flt.xc =repmat((0.5+(0:(flt.N-1)))'/flt.N,1,2).*repmat((B-A)',flt.N,1);

%fault dip and normal vectors 
flt.dv=-[cosd(flt.dip), sind(flt.dip)];
flt.nv= [sind(flt.dip),-cosd(flt.dip)];

% effective confining pressure on fault (MPa)
flt.sigmab=50;

% frictional parameters ( velocity-weakening friction, a-b < 0 )
flt.a=1e-2*ones(flt.N,1);
flt.b=flt.a-5e-3;

% velocity-weakening
vw = abs(flt.xc(:,2))>=.2e3 & abs(flt.xc(:,2))<=35e3;
flt.b(vw) = flt.a(vw) + 5e-3;

% static friction coefficient
flt.mu0=0.6*ones(flt.N,1);

% plate velocity (m/s)
flt.Vpl=1e-9*ones(flt.N,1);

% radiation damping coefficient
flt.damping=5;

% reference slip rate (m/s)
flt.Vo=1e-6*ones(flt.N,1);

%% Shear zone geometry
% consistss of two parts: Oceanic mantle (ocean) and continental mantle (cmantle) 
%% oceanic mantle : when the shear zone is below the depth of megathrust
shzWidth=500e3;
shzdip=flt.dip(1)+2;
Tpl=20e3;
C=[B(1)-shzWidth;B(2)+10e3]; % orgbelow the fault 
D=[ C(1)+shzWidth;   C(2)];% orgbelow the fault 
E=[ D(1)+cosd(shzdip)*shzWidth;  D(2)+sind(shzdip)*shzWidth];
F=[C(1);  E(2)];

% % oceanic mantle : when the shear zone is below Tpl
% shzWidth=500e3;
% Tpl=20e3;
% C=[B(1)-shzWidth;Tpl]; %  
% D=[ A(1);   Tpl];%  
% E=[ D(1)+cosd(flt.dip(1))*shzWidth;  B(2)+10e3+sind(flt.dip(1))*shzWidth];
% F=[C(1);  E(2)];

% GEOMETRY: triangular mesh shear zones: using unicycle distmesh function 
pv=[C,D,E,F,C]';

ocean.meshSize=20e3; %mesh size for shear zones 

figure(1);clf;set(gcf,'name','Mesh optimization (distmesh)')
[ocean.p,ocean.t]=distmesh.distmesh2d(@distmesh.dpoly,@distmesh.huniform,ocean.meshSize, ...
    [min(pv(:,1)),min(pv(:,2)); max(pv(:,1)),max(pv(:,2))],pv,pv);

ocean.N=size(ocean.t,1);

ocean.xc=[(ocean.p(ocean.t(:,1),1)+ocean.p(ocean.t(:,2),1)+ocean.p(ocean.t(:,3),1))/3, ...
          (ocean.p(ocean.t(:,1),2)+ocean.p(ocean.t(:,2),2)+ocean.p(ocean.t(:,3),2))/3];

close(1)
%%
% continental mantle
G=B;
H=[ G(1)+shzWidth;   G(2)];
J=[H(1);  B(2)+10e3+sind(flt.dip(1))*shzWidth];

% GEOMETRY: triangular mesh shear zones: using unicycle distmesh function
pv=[G,H,J,G]';

cmantle.meshSize=20e3;

figure(1);clf;set(gcf,'name','Mesh optimization (distmesh)')
[cmantle.p,cmantle.t]=distmesh.distmesh2d(@distmesh.dpoly,@distmesh.huniform,cmantle.meshSize, ...
    [min(pv(:,1)),min(pv(:,2)); max(pv(:,1)),max(pv(:,2))],pv,pv);

cmantle.N=size(cmantle.t,1);

cmantle.xc=[(cmantle.p(cmantle.t(:,1),1)+cmantle.p(cmantle.t(:,2),1)+cmantle.p(cmantle.t(:,3),1))/3, ...
          (cmantle.p(cmantle.t(:,1),2)+cmantle.p(cmantle.t(:,2),2)+cmantle.p(cmantle.t(:,3),2))/3];

close(1);
%% Rheology of Shear Zones: 1) Oceanic mantle 2) Continental mantle
% All input rheological parameters go here 
% Choose rheology using these options

power = 1;% strain rate = A*stress^(power)
burger = 1;% on-1/off-0
% rheological coefficient
etaval = 1e20;%Maxwell viscosity in Pa-s; for power>1, this is A^{-1}
alphaval = 1/100;% for Kelvin element (eta_K = alpha*etaM)
if burger==1
    power = 1;
    etaM = etaval;
    etaK = etaval*alphaval;
end

% rheological properties ?
% ocean.A = 1./(ocean.A.^(ocean.power)).*Edot.^((1/ocean.power)-1); 

ocean.power=1;
ocean.eta = etaM.*1e-6.*ones(ocean.N,1);
cmantle.eta = etaM.*1e-6.*ones(cmantle.N,1);


% background strain rate tensor (1 / s) gaussian 
ocean_x_epl=linspace(min(ocean.xc(:,1)/1e3),max(ocean.xc(:,1)/1e3),ocean.N);
epl = gaussmf(ocean_x_epl,[40,-100]); % change accordingly
ocean.e22pl= interp1(ocean_x_epl,epl,ocean.xc(:,1)/1e3)*(-1e-15);
ocean.e23pl=zeros(ocean.N,1);
ocean.e33pl=zeros(ocean.N,1);

cmantle_x_epl=linspace(min(cmantle.xc(:,1)/1e3),max(cmantle.xc(:,1)/1e3),cmantle.N);
epl = gaussmf(cmantle_x_epl,[40,450]);
cmantle.e22pl= interp1(cmantle_x_epl,epl,cmantle.xc(:,1)/1e3)*(-1e-15);
cmantle.e23pl=zeros(cmantle.N,1);
cmantle.e33pl=zeros(cmantle.N,1);

%% combining oceanic and continental mantle into one structure here 
mantle.N=ocean.N+cmantle.N;
mantle.p=[ocean.p;cmantle.p];
mantle.t=[ocean.t;cmantle.t+ones(size(cmantle.t))*size(ocean.p,1)];
mantle.xc=[ocean.xc;cmantle.xc];
mantle.e22pl=[ocean.e22pl;cmantle.e22pl];
mantle.e23pl=[ocean.e23pl;cmantle.e23pl];
mantle.e33pl=[ocean.e33pl;cmantle.e33pl];

%% observation points
obs.x=([1;0]*(linspace(-300,750,200)))'*1e3;
obs.N=size(obs.x,1);
obs.dgf=2;  % (east displacement, down displacement in 2D)

%% Geomtery and properties check 
figure(2);clf;set(gcf,'name','Geometry');hold on
plot(flt.x(:,1)/1e3,flt.x(:,2)/1e3)
plot(obs.x(:,1)/1e3,obs.x(:,2)/1e3,'^')
toplot=(mantle.e22pl);
patch([mantle.p(mantle.t(:,1),1),mantle.p(mantle.t(:,2),1),mantle.p(mantle.t(:,3),1),mantle.p(mantle.t(:,1),1)]'/1e3, ...
      [mantle.p(mantle.t(:,1),2),mantle.p(mantle.t(:,2),2),mantle.p(mantle.t(:,3),2),mantle.p(mantle.t(:,1),2)]'/1e3, ...
      [toplot,toplot,toplot,toplot]');
% plot(mantle.xc(mantle.t(1:2,1),1)/1e3,mantle.xc(mantle.t(1:2,1),2)/1e3,'r^')
colorbar
box on
axis equal tight
set(gca,'YDir','Reverse');
xlabel('Distance (km)');
ylabel('Depth (km)');
grid on; 

% return
%% Temp function 
% function erf(x){p=0.3275911;a1=0.254829592;a2=-0.284496736;a3=1.421413741;a4=-1.453152027;a5=1.061405429;t=1/(1+p*x);
% return 1-(a1*t+a2*t**2+a3*t**3+a4*t**4+a5*t**5)*exp(-x**2)}

%% Green's functions for fault elements on megathrust

% Self interactions KK is the shear traction on fault due to fault slip. 
flt.KK=zeros(flt.N);

for k=1:flt.N
 
    T=flt.W(k)/flt.N/1e3;
    [s22,s23,s33]=unicycle.greens.computeStressPlaneStrainShearZone( ...
        flt.xc(:,1),flt.xc(:,2),flt.x(k,1),flt.x(k,2), ...
        T,flt.W(k)/flt.N,flt.dip(k),0,-1/2/T,0,rig,nu);

    t=[s22.*flt.nv(:,1)+s23.*flt.nv(:,2), ...
       s23.*flt.nv(:,1)+s33.*flt.nv(:,2) ];
    ts=t(:,1).*flt.dv(:,1)+t(:,2).*flt.dv(:,2);
    
    flt.KK(:,k)=ts(:);
end
toc

% KL{i} is the stress in the i direction due to fault slip, 
% with i in 22, 23, 33 (in this order).
flt.KL=cell(3,1);

% initialize kernels
for i=1:3
    flt.KL{i}=zeros(mantle.N,flt.N);
end

tic
for k=1:flt.N
    
    
    T=flt.W(k)/flt.N/1e3;
    [s22,s23,s33]=unicycle.greens.computeStressPlaneStrainShearZone( ...
        mantle.xc(:,1),mantle.xc(:,2),flt.x(k,1),flt.x(k,2), ...
        T,flt.W(k)/flt.N,flt.dip(k),0,-1/2/T,0,rig,nu);
        
    flt.KL{1}(:,k)=s22(:);
    flt.KL{2}(:,k)=s23(:);
    flt.KL{3}(:,k)=s33(:);

end
toc

%% KO{i} is the displacement in the i direction due to fault slip, 
% with i (2, 3) (in this order).
flt.KO=cell(obs.dgf,1);

% initialize kernels
for i=1:obs.dgf
    flt.KO{i}=zeros(obs.N,flt.N);
end

tic
for k=1:flt.N
    T=flt.W(k)/flt.N/1e3;
    [u2,u3]=unicycle.greens.computeDisplacementPlaneStrainShearZone( ...
        obs.x(:,1),obs.x(:,2),flt.x(k,1),flt.x(k,2), ...
        T,flt.W(k)/flt.N,flt.dip(k),0,-1/2/T,0,rig,nu);
        
    flt.KO{1}(:,k)=u2(:);
    flt.KO{2}(:,k)=u3(:);
end
toc
%% Green's functions for triangular elements in mantle (using greens function from the bssa 2018)

% LL{i,j} is the stress in the i direction due to strain in the j
% direction, with i,j in 22, 23, 33 (in this order).
mantle.LL=cell(3,3);

% initialize stress kernels
for i=1:3
    for j=1:3
        mantle.LL{i,j}=zeros(mantle.N);
    end
end

% source strain 100;010;001
I=eye(3);

tic
for k=1:mantle.N
    M=mantle.p(mantle.t(k,1),:);
    N=mantle.p(mantle.t(k,2),:);
    P=mantle.p(mantle.t(k,3),:);
    
    for i=1:3
        [s22,s23,s33]=computeStressPlaneStrainTriangleShearZoneFiniteDifference( ...
            mantle.xc(:,1),mantle.xc(:,2),M,N,P,I(i,1),I(i,2),I(i,3),rig,nu);
        
        mantle.LL{1,i}(:,k)=s22(:);
        mantle.LL{2,i}(:,k)=s23(:);
        mantle.LL{3,i}(:,k)=s33(:);
    end
end
fprintf('\n');
toc

% LK{i} is the shear traction on fault due to strain in the i
% direction, with i in 22, 23, 33 (in this order).
mantle.LK=cell(3,1);

% initialize traction kernels
for i=1:3
    mantle.LK{i}=zeros(flt.N,mantle.N);
end

% source strain 100;010;001
I=eye(3);

tic
for k=1:mantle.N
    if 0==mod(fix((k-1)/mantle.N*1000),50)
        fprintf('.');
    end
    M=mantle.p(mantle.t(k,1),:);
    N=mantle.p(mantle.t(k,2),:);
    P=mantle.p(mantle.t(k,3),:);
    
    for i=1:3
        [s22,s23,s33]=computeStressPlaneStrainTriangleShearZoneFiniteDifference( ...
            flt.xc(:,1),flt.xc(:,2),M,N,P,I(i,1),I(i,2),I(i,3),rig,nu);
        
        t=[s22.*flt.nv(:,1)+s23.*flt.nv(:,2), ...
           s23.*flt.nv(:,1)+s33.*flt.nv(:,2) ];
        ts=t(:,1).*flt.dv(:,1)+t(:,2).*flt.dv(:,2);
        
        mantle.LK{i}(:,k)=ts(:);
    end
end
toc

%% LO{i,j} is the displacement component i (i=2,3) at observation points 
% due to strain in the j direction, with j in 22, 23, 33 (in this order).

ocean.LO=cell(obs.dgf,3);

% initialize displacement kernels
for i=1:obs.dgf
    for j=1:3
        ocean.LO{i,j}=zeros(obs.N,ocean.N);
    end
end

% source strain 100;010;001
I=eye(3);

tic
for k=1:ocean.N
    M=ocean.p(ocean.t(k,1),:);
    N=ocean.p(ocean.t(k,2),:);
    P=ocean.p(ocean.t(k,3),:);
    
    
    for i=1:3
        [u2,u3]=computeDisplacementPlaneStrainTriangleShearZone( ...
            obs.x(:,1),obs.x(:,2),M,N,P,I(i,1),I(i,2),I(i,3),nu);
        
        ocean.LO{1,i}(:,k)=u2(:);
        ocean.LO{2,i}(:,k)=u3(:);
    end
end
toc

cmantle.LO=cell(obs.dgf,3);

% initialize displacement kernels
for i=1:obs.dgf
    for j=1:3
        cmantle.LO{i,j}=zeros(obs.N,cmantle.N);
    end
end

% source strain 100;010;001
I=eye(3);

tic
for k=1:cmantle.N
    if 0==mod(fix((k-1)/cmantle.N*1000),50)
        fprintf('.');
    end
    M=cmantle.p(cmantle.t(k,1),:);
    N=cmantle.p(cmantle.t(k,2),:);
    P=cmantle.p(cmantle.t(k,3),:);
    
    
    for i=1:3
        [u2,u3]=computeDisplacementPlaneStrainTriangleShearZone( ...
            obs.x(:,1),obs.x(:,2),M,N,P,I(i,1),I(i,2),I(i,3),nu);
        
        cmantle.LO{1,i}(:,k)=u2(:);
        cmantle.LO{2,i}(:,k)=u3(:);
    end
end
toc
% return
%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                       %
%         N U M E R I C A L   S O L U T I O N           %
%                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
Teq = 200*3.15e7;% earthquake every Teq years
taumax = 3;
eqslip = compute_earthquakeslip_modified(flt,Teq*flt.Vpl(1),taumax);
taueq = flt.KK*eqslip;
taueq(taueq<0) = 0;

if true
    figure(10),clf
    subplot(211)
    plot(flt.xc(:,1)./1e3,eqslip,'LineWidth',2)
    axis tight
    subplot(212)
    plot(flt.xc(:,1)./1e3,flt.KK*eqslip,'rx'), hold on
%     plot(shz.xc(:,1)./1e3,evl.KL*eqslip,'kx')
    axis tight, grid on
    ylim([-0.1 1].*taumax)
end



%%
tic
% state parameters
flt.dgf=2;
mantle.dgf=2;

% initialize State Vector: faults 
Y0=zeros(flt.N*flt.dgf+mantle.N*mantle.dgf,1);
Y0(1:flt.dgf:flt.N*flt.dgf) = zeros(flt.N,1); % slip
Y0(2:flt.dgf:flt.N*flt.dgf) = log(flt.Vpl*0.99./flt.Vo); % log of velocity

% initialize State Vector :shear zones
% stresses -> Refer to Rishav's document 
LLedit = [0.5*(mantle.LL{1,1}-mantle.LL{1,3}-mantle.LL{3,1}+mantle.LL{3,3}) 0.5*(mantle.LL{1,2}-mantle.LL{3,2}); mantle.LL{2,1}-mantle.LL{2,3} mantle.LL{2,2}]; %change to 2n x 2n matrix. remove cells 
% construct an approximation by removing eigen values 
% construct 2x2 matrix based on the right hand side of the equation; 
% 9 components to 4
% eigen values; check and replace greater than 0; 
[V,D] = eig(LLedit);

% change the eigen value to 0 (for positive eig values) to avoid blow ups 
% reconstruct LL.edit T lamba T' 
eig_c = diag(D);
eig_c(real(eig_c)>0) = 0;

bigLmod = V*diag(eig_c)/V;
%  recast as cells to be used in ODE. 
mantle.LLdev=mat2cell(bigLmod,[mantle.N mantle.N],[mantle.N mantle.N]);
%mantle. LL%
% Note here s22 is not the whole stress; its only the deviatoric part;
Y0(flt.N*flt.dgf+1:mantle.dgf:end) = zeros(mantle.N,1); % s22_dev
Y0(flt.N*flt.dgf+2:mantle.dgf:end) = zeros(mantle.N,1); % s23


%% notes
% sigma22_plus = sigma22_minus + KL_22*s 
% sigma22_plus = sigma_22_plus - p_plus

%% Numerical solution 
% initialize function handle with constitutive parameters
yp=@(t,y) imposed_viscoplastic_ode(t,y,flt,mantle);
ncycles=100;
close all

tic
% solve
options=odeset('Refine',1,'AbsTol',1e-6,'RelTol',1e-6,'InitialStep',1e-6,'MaxStep',3e8,'oDir','ode_out/'); 
for i = 1:ncycles
    if i==1
        [sim.t,sim.Y]=ode45(yp,[0 Teq],Y0,options);
    else
        % provide new initial conditions (loaded by earthquake)
        Y0 = sim.Y(end,:);
        
        Y0(1:flt.dgf:flt.N*flt.dgf) = eqslip; 
        Y0(2:flt.dgf:flt.N*flt.dgf) = Y0(2:flt.dgf:flt.N*flt.dgf)' + taueq./(flt.a-flt.b)./flt.sigmab;
        
        % s22dev, s23
        Y0(flt.N*flt.dgf+1:mantle.dgf:end) = Y0(rcv.N*flt.dgf+1:mantle.dgf:end)' + flt.KL{1}-((flt.KL{1}+flt.KL{3})/2)*eqslip; % second term: removing pressure because pressure change (elastic part) from fault slip does not drive nay viscous flow.  
        Y0(flt.N*flt.dgf+2:mantle.dgf:end) = Y0(rcv.N*flt.dgf+2:mantle.dgf:end)' + flt.KL{2}*eqslip; 

        [tmod,Ymod]=ode45(yp,[0,Teq],Y0,options);
        sim.t = tmod;
        sim.Y = Ymod;        
    end
    disp(['Cycle ' num2str(i) ' in progress'])
end

toc

beep

