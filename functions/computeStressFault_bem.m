function [Sxx,Sxz,Szz] = computeStressFault_bem(params,boundary,x,z,nu,G)
% INPUTS
% params = [src.x(i,2) src.x(i,1) src.W(i) src.dip(i) 1];
% boundary - data structure/object containing boundary mesh must contain 
%           xc (patch centers), 
%           nv (unit normal vector for each patch), 
%           sv (unit vector in dip direction)
% x,z - observation points
% nu,G - elastic parameters
% 
% OUTPUTS
% Sxx,Sxz,Szz - stress components after accounting for boundary conditions

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