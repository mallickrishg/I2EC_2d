function [Gx,Gz] = computeFaultDisplacementKernelsBem(src,obs,boundary,scalar_value)
% traction or shear kernel computation for a given source and receiver object pair
% INPUTS
% src - object or data structure containing fault mesh 
%           (end points and center nodes)
% obs - observation points as [x,z] coordinates
% boundary - object or data structure containing boundary mesh for BEM
% 
% OUTPUTS 
% returns half-space kernels 
% ONLY resulting from shear slip on a source  
% 2 x [N x N] matrices containing displacement kernels
% 
% Author:
% Rishav Mallick, JPL, 2023

Nobs = length(obs(:,1));

Gx = zeros(Nobs,src.N);
Gz = zeros(Nobs,src.N);

% compute kernels relating unit shear slip from src to boundary
[Kdd_src_boundary,Kdn_src_boundary,~,~] = geometry.computeFaultTractionKernels(src,boundary);
[Gdx_src_boundary,Gdz_src_boundary,~,~] = geometry.computeFaultDisplacementKernels(src,boundary.xc);

% compute displacement & traction kernels for boundary on itself
[Kdd,Kdn,Knd,Knn] = geometry.computeFaultTractionKernels(boundary,boundary);
[Gdx,Gdz,Gnx,Gnz] = geometry.computeFaultDisplacementKernels(boundary,boundary.xc);
Ktraction = [Kdd,Knd;Kdn,Knn];
Kdisp = [Gdx,Gnx;Gdz,Gnz];

% src to rcv kernels
% and boundary to rcv kernels
[Gdx_src_rcv,Gdz_src_rcv,~,~] = geometry.computeFaultDisplacementKernels(src,obs);
[Gdx_boundary_rcv,Gdz_boundary_rcv,Gnx_boundary_rcv,Gnz_boundary_rcv] = ...
    geometry.computeFaultDisplacementKernels(boundary,obs);

% scalar_value = mean(rcv.Vpl);

for i = 1:src.N
    % we only consider shear sources
    source = zeros(src.N,1);
    source(i) = scalar_value;

    % calculate source tractions
    td_source = Kdd_src_boundary*source;
    tn_source = Kdn_src_boundary*source;
    % calculate source displacements
    ux_source = Gdx_src_boundary*source;
    uz_source = Gdz_src_boundary*source;

    % boundary condition vector
    BC = zeros(2*boundary.N,1);
    % boundary conditions (only 2 types)
    % rcv.Vpl = 0: dirichlet BC (u = 0)
    % rcv.Vpl = 1: Neumann BC (t = 0)
    idirichlet = boundary.Vpl == 0;
    ineumann = boundary.Vpl == 1;

    % displacement BC
    BC([idirichlet;idirichlet]) = 0.*[boundary.Vx(idirichlet);boundary.Vz(idirichlet)] - ...
                                  [ux_source(idirichlet);uz_source(idirichlet)];
    % traction BC
    BC([ineumann;ineumann]) = 0 - [td_source(ineumann);tn_source(ineumann)];        
    
    % construct BEM kernel
    kernel = zeros(boundary.N*2,boundary.N*2);
    kernel([idirichlet;idirichlet],:) = Kdisp([idirichlet;idirichlet],:);
    kernel([ineumann;ineumann],:) = Ktraction([ineumann;ineumann],:);
    
    % BEM solve
    bem_sol = kernel\BC;
    
    % compute displacements from bem solution
    disp = [Gdx_src_rcv;Gdz_src_rcv]*source + ...
        [Gdx_boundary_rcv,Gnx_boundary_rcv;Gdz_boundary_rcv,Gnz_boundary_rcv]*bem_sol;

    Gx(:,i) = disp(1:end/2);
    Gz(:,i) = disp(end/2+1:end);
   
end

end