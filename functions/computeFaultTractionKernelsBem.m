function [varargout] = computeFaultTractionKernelsBem(src,rcv,boundary)
% COMPUTEFAULTTRACTIONKERNELSBEM  Traction or stress kernels due to fault slip, BEM-corrected.
%
%   [Ktau, Ksigma] = computeFaultTractionKernelsBem(src, rcv, boundary)
%   [Kxx, Kzz, Kxz] = computeFaultTractionKernelsBem(src, rcv, boundary)
%
%   For unit dip slip on each source patch, computes the half-space
%   traction (two outputs) or stress (three outputs) at the receiver
%   centres and adds the response to the model boundary, found by a
%   boundary element solve. Only shear (dip-slip) sources are considered.
%
% INPUTS
%   src       geometry.receiver object: source fault patches (Ns patches)
%   rcv       receivers (Nr): a geometry.receiver object for tractions;
%             any object with centres .xc [Nr x 2] (m) for stresses
%   boundary  geometry.receiver object: model boundary patches
%             .Vpl  0: zero displacement; 1: zero traction (.Vx, .Vz
%             are ignored, see NOTES)
%
% OUTPUTS
%   Ktau            [Nr x Ns] shear traction in the receiver dip direction
%                   (MPa/m)
%   Ksigma          [Nr x Ns] normal traction on the receiver (MPa/m)
%   Kxx, Kzz, Kxz   [Nr x Ns] stress components at the receiver centres
%                   (MPa/m)
%
% NOTES
%   - Boundary correction: for each source, slip on the boundary patches
%     is solved so that the total displacement (source plus boundary) is
%     zero where boundary.Vpl == 0 (Dirichlet) and the total traction is
%     zero where boundary.Vpl == 1 (Neumann).
%   - The kernels solve the homogeneous problem only. The prescribed
%     boundary motion (.Vx, .Vz, from bcx, bcz in the .seg file) is
%     ignored; it would be a separate particular solution added once, not
%     part of each column. The background rate comes from the long-term
%     solution instead.
%   - Stress components are returned in the order (xx, zz, xz), unlike
%     computeShzStressKernelsBem, which uses (xx, xz, zz).
%   - Errors unless called with 2 or 3 outputs.
%
% SEE ALSO
%   computeAllStressKernelsBem, computeFaultDisplacementKernelsBem,
%   computeShzStressKernelsBem
%
% AUTHORS
%   Rishav Mallick (JPL), 2023

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
[Kdd_src_boundary,Kdn_src_boundary,~,~] = geometry.computeFaultTractionKernels(src,boundary);
[Gdx_src_boundary,Gdz_src_boundary,~,~] = geometry.computeFaultDisplacementKernels(src,boundary.xc);

% compute displacement & traction kernels for boundary on itself
[Kdd,Kdn,Knd,Knn] = geometry.computeFaultTractionKernels(boundary,boundary);
[Gdx,Gdz,Gnx,Gnz] = geometry.computeFaultDisplacementKernels(boundary,boundary.xc);
Ktraction = [Kdd,Knd;Kdn,Knn];
Kdisp = [Gdx,Gnx;Gdz,Gnz];

% src to rcv kernels
% and boundary to rcv kernels
if nargout == 2
    [Kdd_src_rcv,Kdn_src_rcv,~,~] = geometry.computeFaultTractionKernels(src,rcv);
    [Kdd_boundary_rcv,Kdn_boundary_rcv,Knd_boundary_rcv,Knn_boundary_rcv] = ...
        geometry.computeFaultTractionKernels(boundary,rcv);
elseif nargout == 3
    [Kdxx_src_rcv,Kdzz_src_rcv,Kdxz_src_rcv,~,~,~] = geometry.computeFaultStressKernels(src,rcv);
    [Kdxx_boundary_rcv,Kdzz_boundary_rcv,Kdxz_boundary_rcv,Knxx_boundary_rcv,Knzz_boundary_rcv,Knxz_boundary_rcv] = geometry.computeFaultStressKernels(boundary,rcv);
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

    % displacement BC: homogeneous (boundary.Vx, Vz are not used)
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

%  outputs
if nargout == 2 % fault traction kernels
    varargout{1} = Ktau;
    varargout{2} = Ksigma;
elseif nargout == 3 % full stress kernels
    varargout{1} = Kxx;
    varargout{2} = Kzz;
    varargout{3} = Kxz;
else
end

end