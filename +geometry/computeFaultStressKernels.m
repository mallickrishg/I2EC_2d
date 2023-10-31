function [Kdxx,Kdzz,Kdxz,Knxx,Knzz,Knxz] = computeFaultStressKernels(src,rcv)
% full stress kernel computation in 2-d plane strain
% returns all 6 stress kernels where
% each kernel is a component of the full 2-d stress tensor [sxx,szz,sxz]
% 
% OUTPUTS:
% Kdxx,Kdzz,Kdxz - for dip-slip source component
% Knxx,Knzz,Knxz - for tensile-slip source component
% 
% Author:
% Rishav Mallick, JPL, 2023

Kdxx = zeros(rcv.N,src.N);
Kdzz = zeros(rcv.N,src.N);
Kdxz = zeros(rcv.N,src.N);
Knxx = zeros(rcv.N,src.N);
Knzz = zeros(rcv.N,src.N);
Knxz = zeros(rcv.N,src.N);

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
    
    Kdxx(:,i) = Sxx;
    Kdzz(:,i) = Szz;
    Kdxz(:,i) = -Sxz;
    
    % for tensile component
    [Stress] = geometry.LDstressHS(x,z,src.xc(i,1),src.xc(i,2),src.W(i)/2,-deg2rad(src.dip(i)),...
        0,1,src.earthModel.nu,2*src.earthModel.G*(1+src.earthModel.nu));
    Sxx = Stress(:,1);
    Szz = Stress(:,2);
    Sxz = Stress(:,3);
    
    Knxx(:,i) = Sxx;
    Knzz(:,i) = Szz;
    Knxz(:,i) = -Sxz;
end



end