function [varargout] = computeFaultTractionKernelsBem(src,rcv,boundary)
% traction kernel computation for a given source and receiver object pair
% INPUTS
% src,rcv - objects or data structures containing fault mesh 
%           (end points and center nodes)
% boundary - object or data structure containing boundary mesh for BEM
% 
% OUTPUTS 
% returns half-space kernels 
% ONLY resulting from shear slip on a source 
% based on number of output variables
% 
% if varargout == 2 (traction kernels)
%   2 x [N x N] matrices containing fault-centric traction kernels
%   Ktau - traction kernel in shear direction
%   Ksigma - traction kernel in fault-normal direction
% if varargout == 3 (full stress kernel)
%   3 x [N x N] matrices containing stress kernels
%   Kxx,Kzz,Kxz
% 
% Author:
% Rishav Mallick, JPL, 2023

if nargout == 2
    Ktau = zeros(rcv.N,src.N);
    Ksigma = zeros(rcv.N,src.N);
elseif nargout == 3
    Kxx = zeros(rcv.N,src.N);
    Kzz = zeros(rcv.N,src.N);
    Kxz = zeros(rcv.N,src.N);
else
    error('incorrect number of output options. provide either 2 (traction) or 3 (stress) output variables.')
end

% compute kernels relating unit shear slip from src to boundary
[Kdd_src_boundary,Kdn_src_boundary,~,~] = geometry.computeFullTractionKernels(src,boundary);
[Gdx_src_boundary,Gdz_src_boundary,~,~] = geometry.computeDisplacementKernels(src,boundary.xc);

% compute displacement & traction kernels for boundary on itself
[Kdd,Kdn,Knd,Knn] = geometry.computeFullTractionKernels(boundary,boundary);
[Gdx,Gdz,Gnx,Gnz] = geometry.computeDisplacementKernels(boundary,boundary.xc);
Ktraction = [Kdd,Knd;Kdn,Knn];
Kdisp = [Gdx,Gnx;Gdz,Gnz];

% src to rcv kernels
% and boundary to rcv kernels
if nargout == 2
    [Kdd_src_rcv,Kdn_src_rcv,~,~] = geometry.computeFullTractionKernels(src,rcv);
    [Kdd_boundary_rcv,Kdn_boundary_rcv,Knd_boundary_rcv,Knn_boundary_rcv] = ...
        geometry.computeFullTractionKernels(boundary,rcv);
elseif nargout == 3
    [Kdxx_src_rcv,Kdzz_src_rcv,Kdxz_src_rcv,~,~,~] = geometry.computeFullStressKernels(src,rcv);
    [Kdxx_boundary_rcv,Kdzz_boundary_rcv,Kdxz_boundary_rcv,Knxx_boundary_rcv,Knzz_boundary_rcv,Knxz_boundary_rcv] = geometry.computeFullStressKernels(boundary,rcv);
else
end

% boundary to rcv kernels

for i = 1:src.N
    % we only consider shear sources
    source = zeros(src.N,1);
    source(i) = 1;

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
    BC([idirichlet;idirichlet]) = 0 - [ux_source(idirichlet);uz_source(idirichlet)];
    % traction BC
    BC([ineumann;ineumann]) = 0 - [td_source(ineumann);tn_source(ineumann)];        
    
    % construct BEM kernel
    kernel = zeros(boundary.N*2,boundary.N*2);
    kernel([idirichlet;idirichlet],:) = Kdisp([idirichlet;idirichlet],:);
    kernel([ineumann;ineumann],:) = Ktraction([ineumann;ineumann],:);
    
    % BEM solve
    bem_sol = kernel\BC;
    
    % compute tractions or stress from bem solution
    
    if nargout == 2
        tau = [Kdd_src_rcv;Kdn_src_rcv]*source + ...
            [Kdd_boundary_rcv,Knd_boundary_rcv;Kdn_boundary_rcv,Knn_boundary_rcv]*bem_sol;

        Ktau(:,i) = tau(1:end/2);
        Ksigma(:,i) = tau(end/2+1:end);        
    elseif nargout == 3
        sigma = [Kdxx_src_rcv;Kdzz_src_rcv;Kdxz_src_rcv]*source + ...
            [Kdxx_boundary_rcv,Knxx_boundary_rcv;Kdzz_boundary_rcv,Knzz_boundary_rcv;Kdxz_boundary_rcv,Knxz_boundary_rcv]*bem_sol;
        Kxx(:,i) = sigma(1:end/3);
        Kzz(:,i) = sigma(end/3+1:2/3*end);
        Kxz(:,i) = sigma(2/3*end+1:end);
    else

    end
end

if nargout == 2
    varargout{1} = Ktau;
    varargout{2} = Ksigma;
elseif nargout == 3
    varargout{1} = Kxx;
    varargout{2} = Kzz;
    varargout{3} = Kxz;
else
end

end