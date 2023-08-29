% script to create a triangular mesh given N (x,z) points
% for oceanic and continental mantle
% Rishav Mallick, Caltech Seismolab, 2023

clear

% mesh close to a given point
% specify location (x0,y0)
x0 = 2000;
y0 = -500;
% create a function that describes how the mesh should vary from x0,y0
hmin = 1;
dh = 2;

hfun1 = @(x,y) hmin + dh.*sqrt((x - x0).^2  + (y - y0).^2);
%% create oceanic mantle mesh

% provide 5 points to 
pt1 = [-100,-20];
pt2 = [0,-20];
pt3 = [200,-70];
pt4 = [200,-80];
pt5 = [-100,-80];

v_o = [pt1;...
    pt2;...
    pt3;...
    pt4;...
    pt5];

hdata = [];
hdata.hmax = 10;
hdata.fun = hfun1;

% this is the meshing routine
[p_o, t_o] = mesh2d ( v_o, [], hdata);
% it creates a figure by default
close


figure(1),clf
plot(v_o(:,1),v_o(:,2),'ko','MarkerFaceColor','r')
hold on
plot([0,200],[0,-50],'ko-','LineWidth',3)
for i = 1:length(t_o)
    plot([p_o(t_o(i,:),1);p_o(t_o(i,1),1)],[p_o(t_o(i,:),2);p_o(t_o(i,1),2)],'r-')
end
axis equal
set(gca,'YDir','normal','FOntsize',20,'Linewidth',2)

%% create second mesh for continental mantle

% provide 5 points to 
pt1 = [201,-50];
pt2 = [201,-35];
pt3 = [350,-35];
pt4 = [350,-75];
pt5 = [300,-75];

v_c = [pt1;...
    pt2;...
    pt3;...
    pt4;...
    pt5];

% this is the meshing routine
[p_c, t_c] = mesh2d ( v_c, [], hdata);
% it creates a figure by default
close

figure(1)
plot(v_c(:,1),v_c(:,2),'ko','MarkerFaceColor','r')
hold on
for i = 1:length(t_c)
    plot([p_c(t_c(i,:),1);p_c(t_c(i,1),1)],[p_c(t_c(i,:),2);p_c(t_c(i,1),2)],'r-')
end
axis equal
set(gca,'YDir','normal','Fontsize',20,'Linewidth',2)

%% save both meshes in a text file
p = [p_o;p_c];
t = [t_o;t_c + length(p_o(:,1))];

figure(2),clf
plot([0,200],[0,-50],'ko-','LineWidth',3)
hold on
for i = 1:length(t)
    plot([p(t(i,:),1);p(t(i,1),1)],[p(t(i,:),2);p(t(i,1),2)],'r-')
end
axis tight equal
set(gca,'YDir','normal','FOntsize',20,'Linewidth',2)

writetable(table(p.*1e3),'../shearzone_vertices.dat','WriteVariableNames',false)
writetable(table(t),'../shearzone_triangulation.dat','WriteVariableNames',false)





