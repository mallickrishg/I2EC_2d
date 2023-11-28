clear
% close all

% mesh close to a given point
% specify location (x0,y0)
x0 = 40;
y0 = 50;
% create a function that describes how the mesh should vary from x0,y0
hmin = 20;
dh = 0;

hfun1 = @(x,y) hmin + dh*sqrt(( x - x0 ).^2  + ( y - y0 ).^2);

% v =[0.0, 0.0;...
%     2.0, 0.0;...
%     2.0, 1.0;...
%     2.0, 2.0;...
%     1.0, 2.0;...
%     0.0, 2.0 ];

v = [0,20;...
    210,20;...
    210,100;...
    0,50];

hdata = [];
hdata.fun = hfun1;

% this is the meshing routine
[ p, t ] = mesh2d ( v, [], hdata);
% it creates a figure by default
close


figure(1),clf
plot(v(:,1),v(:,2),'ko','MarkerFaceColor','r')
hold on
for i = 1:length(t)
    plot([p(t(i,:),1);p(t(i,1),1)],[p(t(i,:),2);p(t(i,1),2)],'r-')
end
axis tight equal

% use MATLAB's meshing routine
% R1 = [3,4,-1,1,1,-1,0.5,0.5,-0.75,-0.75]';
model = createpde(1);
% R1 = [2,4,0,210,210,0,20,20,100,50]';
R1 = [2,length(v(:,1)),v(:,1)',v(:,2)']';

g = decsg(R1);
geometryFromEdges(model,g);
generateMesh(model,"Hmax",12,"GeometricOrder","linear");
figure(1)
pdeplot(model)
axis tight equal
set(gca,'YDir','reverse','FOntsize',20,'Linewidth',2)
