function [Gdx,Gdz,Gnx,Gnz] = computeFaultDisplacementKernels(rcv,ox)
% COMPUTEFAULTDISPLACEMENTKERNELS  Half-space displacements from fault slip.
%
%   [Gdx, Gdz, Gnx, Gnz] = geometry.computeFaultDisplacementKernels(rcv, ox)
%
%   Displacement at the observation points due to unit dip slip and unit
%   opening on each patch of rcv, in a homogeneous elastic half-space
%   (geometry.LDdispHS). No boundary-element correction is applied; see
%   computeFaultDisplacementKernelsBem for that.
%
% INPUTS
%   rcv     geometry.receiver object: source patches
%           .N           [1 x 1] number of patches
%           .xc          [N x 2] patch centres (m, z negative)
%           .W           [N x 1] patch widths (m)
%           .dip         [N x 1] patch dips (degrees)
%           .nv          [N x 2] unit normal vectors
%           .earthModel  geometry.LDhs object (uses .nu)
%   ox      [Nobs x 2] observation points (x, z) (m, z <= 0)
%
% OUTPUTS
%   Gdx, Gdz  [Nobs x N] x and z displacement per unit dip slip (m/m)
%   Gnx, Gnz  [Nobs x N] x and z displacement per unit opening (m/m)
%
% NOTES
%   - An observation point that coincides with a patch centre is moved
%     by 1e-9 * mean(ox(:,1)) along a patch normal. That normal is
%     rcv.nv indexed by the observation point, so this only works when
%     ox are the centres of rcv itself (e.g. boundary.xc for boundary).
%
% SEE ALSO
%   geometry.LDdispHS, computeFaultDisplacementKernelsBem,
%   geometry.computedisplacementGFs
%
% AUTHORS
%   Rishav Mallick, 2023

delta = 1e-9.*mean(ox(:,1));
xobs = ox(:,1);
zobs = ox(:,2);

Gdx = zeros(length(xobs),rcv.N);
Gdz = zeros(length(xobs),rcv.N);
Gnx = zeros(length(xobs),rcv.N);
Gnz = zeros(length(xobs),rcv.N);

for k = 1:rcv.N
    xobs_mod = xobs;
    zobs_mod = zobs;
    r = sqrt((xobs-rcv.xc(k,1)).^2 + (zobs-rcv.xc(k,2)).^2);
    xobs_mod(r==0) = xobs_mod(r==0) + rcv.nv(r==0,1)*delta;
    zobs_mod(r==0) = zobs_mod(r==0) + rcv.nv(r==0,2)*delta;
    
    [Disp] = geometry.LDdispHS(xobs_mod,zobs_mod,rcv.xc(k,1),rcv.xc(k,2),rcv.W(k)/2,...
        -deg2rad(rcv.dip(k)),1,0,rcv.earthModel.nu);
    
    Gdx(:,k) = Disp(:,1);
    Gdz(:,k) = Disp(:,2);
    
    [Disp] = geometry.LDdispHS(xobs_mod,zobs_mod,rcv.xc(k,1),rcv.xc(k,2),rcv.W(k)/2,...
        -deg2rad(rcv.dip(k)),0,1,rcv.earthModel.nu);
    
    Gnx(:,k) = Disp(:,1);
    Gnz(:,k) = Disp(:,2);
end


end