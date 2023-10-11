function [Kdd,Kdn,Knd,Knn] = computeFaultTractionKernels(src,rcv)
% full traction kernel computation in 2-d plane strain
% returns all 4 traction kernels 
% Kdd,Kdn - for dip-slip coponent
% Knd,Knn - for tensile-slip component
% Rishav Mallick, 2023

Kdd = zeros(rcv.N,src.N);
Kdn = zeros(rcv.N,src.N);
Knd = zeros(rcv.N,src.N);
Knn = zeros(rcv.N,src.N);

% observation points
x = rcv.xc(:,1);
z = rcv.xc(:,2);

parfor i = 1:src.N
    
    % for dip-slip component
    [Stress] = geometry.LDstressHS(x,z,src.xc(i,1),src.xc(i,2),src.W(i)/2,-deg2rad(src.dip(i)),...
        1,0,src.earthModel.nu,2*src.earthModel.G*(1+src.earthModel.nu));
    Sxx = Stress(:,1);
    Szz = Stress(:,2);
    Sxz = Stress(:,3);
    
    t=[Sxx.*rcv.nv(:,1)+Sxz.*rcv.nv(:,2), ...
        Sxz.*rcv.nv(:,1)+Szz.*rcv.nv(:,2)];
    
    Kdd(:,i) = rcv.dv(:,1).*t(:,1) + rcv.dv(:,2).*t(:,2);
    Kdn(:,i) = rcv.nv(:,1).*t(:,1) + rcv.nv(:,2).*t(:,2);
    
    % for tensile component
    [Stress] = geometry.LDstressHS(x,z,src.xc(i,1),src.xc(i,2),src.W(i)/2,-deg2rad(src.dip(i)),...
        0,1,src.earthModel.nu,2*src.earthModel.G*(1+src.earthModel.nu));
    Sxx = Stress(:,1);
    Szz = Stress(:,2);
    Sxz = Stress(:,3);
    
    t=[Sxx.*rcv.nv(:,1)+Sxz.*rcv.nv(:,2), ...
        Sxz.*rcv.nv(:,1)+Szz.*rcv.nv(:,2)];
    
    Knd(:,i) = rcv.dv(:,1).*t(:,1) + rcv.dv(:,2).*t(:,2);
    Knn(:,i) = rcv.nv(:,1).*t(:,1) + rcv.nv(:,2).*t(:,2);
end



end