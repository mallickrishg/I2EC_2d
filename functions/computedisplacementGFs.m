function Gd = computedisplacementGFs(rcv,ox)
% COMPUTEDISPLACEMENTGFS  Half-space displacement Green's functions for fault slip.
%
%   Gd = computedisplacementGFs(rcv, ox)
%
%   For unit dip slip on each patch, computes the displacement at the
%   observation points in a homogeneous half-space with LDdispHS. There is
%   no boundary correction.
%
% INPUTS
%   rcv     geometry.receiver object: source patches (N patches)
%           .xc   [N x 2] patch centres [x, z] (m, z < 0)
%           .W    [N x 1] patch widths (m)
%           .dip  [N x 1] dips (degrees)
%           .earthModel  elastic parameters (.nu is used)
%   ox      [Nobs x 2] observation points [x, z] (m, z <= 0)
%
% OUTPUTS
%   Gd      [2*Nobs x N] displacement per unit slip: odd rows ux, even
%           rows uz
%
% NOTES
%   - Not used by the I2EC solvers (the call in
%     computeAllDisplacementKernelsBem is commented out); see
%     computeFaultDisplacementKernelsBem for the BEM-corrected kernels.
%   - The +geometry package has a copy that calls geometry.LDdispHS.
%
% SEE ALSO
%   LDdispHS, computeFaultDisplacementKernelsBem
%
% AUTHORS
%   Rishav Mallick, 2023

xobs = ox(:,1);
zobs = ox(:,2);

Gd = zeros(length(xobs)*2,rcv.N);

for k = 1:rcv.N
    
    [Disp] = LDdispHS(xobs,zobs,rcv.xc(k,1),rcv.xc(k,2),rcv.W(k)/2,...
        -deg2rad(rcv.dip(k)),1,0,rcv.earthModel.nu);
    
    Gd(1:2:end,k) = Disp(:,1);
    Gd(2:2:end,k) = Disp(:,2);
end


end