% script to create a mesh given 4 (x,z) points
% Rishav Mallick, Caltech Seismolab, 2023

clear

pt1 = [-100,-20];
pt2 = [0,-20];
pt3 = [200,-50];
pt4 = [200,-60];
pt5 = [-100,-60];

% mesh close to a given point
% specify location (x0,y0)
x0 = 2000;
y0 = -50;
% create a function that describes how the mesh should vary from x0,y0
hmin = 0.1;
dh = 1;

hfun1 = @(x,y) hmin + dh.*sqrt((x - x0).^2  + (y - y0).^2);


v = [pt1;...
    pt2;...
    pt3;...
    pt4;...
    pt5];

hdata = [];
hdata.hmax = 6.3;
hdata.fun = hfun1;

% this is the meshing routine
[p, t] = mesh2d ( v, [], hdata);
% it creates a figure by default
close


figure(1),clf
plot(v(:,1),v(:,2),'ko','MarkerFaceColor','r')
hold on
for i = 1:length(t)
    plot([p(t(i,:),1);p(t(i,1),1)],[p(t(i,:),2);p(t(i,1),2)],'r-')
end
axis equal
set(gca,'YDir','normal','FOntsize',20,'Linewidth',2)