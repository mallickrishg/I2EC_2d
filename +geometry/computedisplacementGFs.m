function Gd = computedisplacementGFs(rcv,ox)
% COMPUTEDISPLACEMENTGFS  Half-space displacements from dip slip, interleaved.
%
%   Gd = geometry.computedisplacementGFs(rcv, ox)
%
%   Displacement at the observation points due to unit dip slip on each
%   patch of rcv, in a homogeneous elastic half-space (geometry.LDdispHS),
%   with the x and z components interleaved by row. Not called by the
%   I2EC drivers; computeFaultDisplacementKernels returns the same dip-slip
%   kernels as separate matrices.
%
% INPUTS
%   rcv     geometry.receiver object: source patches
%           .N           [1 x 1] number of patches
%           .xc          [N x 2] patch centres (m, z negative)
%           .W           [N x 1] patch widths (m)
%           .dip         [N x 1] patch dips (degrees)
%           .earthModel  geometry.LDhs object (uses .nu)
%   ox      [Nobs x 2] observation points (x, z) (m, z <= 0)
%
% OUTPUTS
%   Gd      [2*Nobs x N] displacement per unit dip slip (m/m): rows
%           1, 3, 5, ... are x, rows 2, 4, 6, ... are z
%
% SEE ALSO
%   geometry.LDdispHS, geometry.computeFaultDisplacementKernels
%
% AUTHORS
%   Rishav Mallick, 2023

xobs = ox(:,1);
zobs = ox(:,2);

Gd = zeros(length(xobs)*2,rcv.N);

for k = 1:rcv.N
    
    [Disp] = geometry.LDdispHS(xobs,zobs,rcv.xc(k,1),rcv.xc(k,2),rcv.W(k)/2,...
        -deg2rad(rcv.dip(k)),1,0,rcv.earthModel.nu);
    
    Gd(1:2:end,k) = Disp(:,1);
    Gd(2:2:end,k) = Disp(:,2);
end


end