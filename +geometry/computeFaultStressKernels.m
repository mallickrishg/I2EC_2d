function [Kdxx,Kdzz,Kdxz,Knxx,Knzz,Knxz] = computeFaultStressKernels(src,rcv)
% COMPUTEFAULTSTRESSKERNELS  Half-space stress from fault slip and opening.
%
%   [Kdxx, Kdzz, Kdxz, Knxx, Knzz, Knxz] = computeFaultStressKernels(src, rcv)
%
%   Full 2-D stress tensor at the centres of rcv due to unit dip slip and
%   unit opening on each patch of src, in a homogeneous elastic
%   half-space (geometry.LDstressHS). Used for fault-to-mantle stress
%   kernels, so rcv is usually a geometry.shearZoneReceiver.
%
% INPUTS
%   src     geometry.receiver object: source patches
%           .N           [1 x 1] number of patches
%           .xc          [N x 2] patch centres (m, z negative)
%           .W           [N x 1] patch widths (m)
%           .dip         [N x 1] patch dips (degrees)
%           .earthModel  geometry.LDhs object (uses .G, .nu)
%   rcv     object with receiver points (geometry.receiver or
%           geometry.shearZoneReceiver)
%           .N           [1 x 1] number of receivers
%           .xc          [Nr x 2] receiver points (m, z negative)
%
% OUTPUTS
%   Kdxx, Kdzz, Kdxz  [Nr x N] stress sxx, szz, sxz per unit dip slip
%                     (MPa/m)
%   Knxx, Knzz, Knxz  [Nr x N] stress sxx, szz, sxz per unit opening
%                     (MPa/m)
%
% NOTES
%   - The component order is [sxx, szz, sxz]; computeShzStressKernels uses
%     [sxx, sxz, szz].
%   - The sign of the shear stress from geometry.LDstressHS is flipped
%     (Kdxz = -Sxy), but the normal stresses are not.
%   - Young's modulus is passed to geometry.LDstressHS as 2*G*(1 + nu).
%
% SEE ALSO
%   geometry.LDstressHS, geometry.computeFaultTractionKernels,
%   computeFaultTractionKernelsBem
%
% AUTHORS
%   Rishav Mallick (JPL), 2023

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