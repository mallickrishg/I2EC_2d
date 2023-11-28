% script to create a triangular mesh given N (x,z) points
% for oceanic and continental mantle using MATLAB's meshing routines
% Rishav Mallick, Caltech Seismolab, 2023

clear

Hmax1 = 10;
Hmax2 = 6;
%% create oceanic mantle mesh
model1 = createpde;

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

R1 = [2,length(v_o(:,1)),v_o(:,1)',v_o(:,2)']';

g1 = decsg(R1);
geometryFromEdges(model1,g1);
mesh1 = generateMesh(model1,"Hmax",Hmax1,"GeometricOrder","linear");


%% create second mesh for continental mantle
model2 = createpde;
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

R2 = [2,length(v_c(:,1)),v_c(:,1)',v_c(:,2)']';

g2 = decsg(R2);
geometryFromEdges(model2,g2);
mesh2 = generateMesh(model2,"Hmax",Hmax2,"GeometricOrder","linear");

%% save both meshes in a text file
figure(1),clf
pdeplot(mesh1), hold on
pdeplot(mesh2)
axis tight equal