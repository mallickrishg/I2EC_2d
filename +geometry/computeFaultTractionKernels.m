function [Kdd,Kdn,Knd,Knn] = computeFaultTractionKernels(src,rcv)
% COMPUTEFAULTTRACTIONKERNELS  Half-space tractions on faults from fault slip.
%
%   [Kdd, Kdn, Knd, Knn] = geometry.computeFaultTractionKernels(src, rcv)
%
%   Shear and normal traction at the centres of rcv due to unit dip slip
%   and unit opening on each patch of src, in a homogeneous elastic
%   half-space. The stress from geometry.LDstressHS is projected onto the
%   dip and normal vectors of each receiver patch.
%
% INPUTS
%   src     geometry.receiver object: source patches
%           .N           [1 x 1] number of patches
%           .xc          [N x 2] patch centres (m, z negative)
%           .W           [N x 1] patch widths (m)
%           .dip         [N x 1] patch dips (degrees)
%           .earthModel  geometry.LDhs object (uses .G, .nu)
%   rcv     geometry.receiver object: receiver patches
%           .N           [1 x 1] number of patches
%           .xc          [Nr x 2] patch centres (m, z negative)
%           .dv          [Nr x 2] unit dip vectors
%           .nv          [Nr x 2] unit normal vectors
%
% OUTPUTS
%   Kdd     [Nr x N] shear traction per unit dip slip (MPa/m)
%   Kdn     [Nr x N] normal traction per unit dip slip (MPa/m)
%   Knd     [Nr x N] shear traction per unit opening (MPa/m)
%   Knn     [Nr x N] normal traction per unit opening (MPa/m)
%
% NOTES
%   - Young's modulus is passed to geometry.LDstressHS as 2*G*(1 + nu).
%
% SEE ALSO
%   geometry.LDstressHS, geometry.computeFaultStressKernels,
%   computeFaultTractionKernelsBem
%
% AUTHORS
%   Rishav Mallick, 2023

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