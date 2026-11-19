function [Sxx,Sxz,Szz] = computeStressFault_bem(params,boundary,x,z,nu,G)
% COMPUTESTRESSFAULT_BEM  Unfinished: edge-dislocation stress with a traction-free boundary.
%
%   [Sxx, Sxz, Szz] = computeStressFault_bem(params, boundary, x, z, nu, G)
%
%   Meant to compute the stress at (x, z) due to an edge dislocation,
%   corrected for a traction-free model boundary. The function is
%   unfinished and does not run (see NOTES); nothing in I2EC calls it.
%
% INPUTS
%   params    [1 x 5] edge dislocation as in EdgeStress, e.g.
%             [src.x(i,2) src.x(i,1) src.W(i) src.dip(i) 1]
%   boundary  geometry.receiver object: boundary patches
%             .xc  [Nb x 2] patch centres (m)
%             .nv  [Nb x 2] unit normal vectors
%             .dv  [Nb x 2] unit vectors in the dip direction
%   x, z      observation points (m)
%   nu        Poisson's ratio
%   G         shear modulus (MPa)
%
% OUTPUTS
%   Sxx, Sxz, Szz  stress components (never assigned; see NOTES)
%
% NOTES
%   - Calls computeFullTractionKernels, which does not exist in the
%     repository, uses an undefined variable ts, and never assigns its
%     outputs. A TODO in the code says the boundary conditions are not yet
%     implemented.
%   - computeFaultTractionKernelsBem computes BEM-corrected stress from
%     fault slip.
%
% SEE ALSO
%   EdgeStress, computeFaultTractionKernelsBem
%
% AUTHORS
%   Rishav Mallick, 2023

% first construct traction kernels for boundary mesh
[Kdd,Kdn,Knd,Knn] = computeFullTractionKernels(boundary,boundary);

% compute stress change from source onto boundary mesh
[Sxx_source,Sxz_source,Szz_source] = EdgeStress(params,boundary.xc(:,1),boundary.xc(:,2),nu,G);
% compute tractions on fault
t = [Sxx_source.*boundary.nv(:,1) + Sxz_source.*boundary.nv(:,2), ...
     Sxz_source.*boundary.nv(:,1) + Szz_source.*boundary.nv(:,2)];
% rotate tractions to fault-centric coordinates
td = boundary.dv(:,1).*t(:,1) + boundary.dv(:,2).*t(:,2);
tn = boundary.nv(:,1).*t(:,1) + boundary.nv(:,2).*t(:,2);

% Solve BEM problem - assuming the boundary is traction free
kernel = [Kdd, Knd;...
          Kdn, Knn];
bem_slip = -pinv(kernel)*[ts;tn];

% TODO - need to implement appropriate disp,traction boundary conditions 
% and create and arrange kernel accordingly


end