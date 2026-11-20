function LL_bem = computeShzStressKernelsBem(src,shz,boundary,bc)
% COMPUTESHZSTRESSKERNELSBEM  Stress or traction kernels due to shear-zone strain, BEM-corrected.
%
%   LL_bem = computeShzStressKernelsBem(src, shz, boundary)
%   LL_bem = computeShzStressKernelsBem(src, shz, boundary, bc)
%
%   For a unit eigenstrain in each source triangle, one strain component at
%   a time, computes the half-space stress (receiver is a
%   geometry.shearZoneReceiver) or traction (receiver is a
%   geometry.receiver) at the receiver centres and adds the response to the
%   model boundary, found by a boundary element solve.
%
% INPUTS
%   src       geometry.shearZoneReceiver object: source triangles (Ns)
%   shz       receivers (Nr): geometry.shearZoneReceiver object (stress) or
%             geometry.receiver object (traction)
%   boundary  geometry.receiver object: model boundary patches
%             .Vpl  0: zero displacement; 1: zero traction (.Vx, .Vz
%             are ignored, see NOTES)
%   bc        (optional) boundary correction of src from
%             computeShzBoundaryCorrection(src, boundary); computed here
%             if omitted or empty. Pass it when building kernels for
%             several receiver sets with the same sources and boundary
%
% OUTPUTS
%   LL_bem    [Nr x Ns x 3 x 3] kernels (MPa per unit strain); 4th index =
%             source strain component (e22, e23, e33), 3rd index =
%             receiver component: (sxx, sxz, szz) for a shearZoneReceiver;
%             (shear traction, normal traction, zeros) for a receiver
%
% NOTES
%   - Source strains are tensor components: e23 = 1 is a unit tensor shear
%     strain.
%   - The half-space stress comes from geometry.computeShzStressKernels
%     (closed form, computeStressPlaneStrainTriangleShearZone).
%   - Boundary correction as in computeFaultTractionKernelsBem: the
%     boundary slip that cancels each source on the boundary
%     (computeShzBoundaryCorrection) depends only on the sources and the
%     boundary, so only the receiver kernels are computed per call when bc
%     is passed. Homogeneous problem only: the prescribed .Vx, .Vz are
%     ignored (see computeFaultTractionKernelsBem).
%   - For a shearZoneReceiver, the boundary contribution from
%     geometry.computeFaultStressKernels is stacked in the same component
%     order as the half-space part, (sxx, sxz, szz).
%
% SEE ALSO
%   computeShzBoundaryCorrection, computeAllStressKernelsBem,
%   computeFaultTractionKernelsBem, computeShzDisplacementKernelsBem
%
% AUTHORS
%   Rishav Mallick (JPL), 2023


if nargin < 4 || isempty(bc)
    bc = computeShzBoundaryCorrection(src,boundary);
elseif bc.Ns ~= src.N || bc.Nb ~= boundary.N
    error('bc was computed for %d sources and %d boundary patches, not %d and %d.', ...
        bc.Ns, bc.Nb, src.N, boundary.N);
end

% half-space stress or traction from the sources at the receivers
LL_src_shz = geometry.computeShzStressKernels(src,shz);

% boundary to receiver kernels, stacked in the output component order
if isa(shz,'geometry.shearZoneReceiver')
    [Kdxx,Kdzz,Kdxz,Knxx,Knzz,Knxz] = geometry.computeFaultStressKernels(boundary,shz);
    Kb = {[Kdxx,Knxx]; [Kdxz,Knxz]; [Kdzz,Knzz]}; % sxx, sxz, szz
elseif isa(shz,'geometry.receiver')
    [Kdd,Kdn,Knd,Knn] = geometry.computeFaultTractionKernels(boundary,shz);
    Kb = {[Kdd,Knd]; [Kdn,Knn]}; % shear, normal traction
else
    error('not a recognized geometry, provide either geometry.receiver or geometry.shearZoneReceiver')
end

% half-space part plus the boundary slip that cancels each source on the
% boundary (a shearZoneReceiver fills all 3 receiver components, a
% receiver the first 2)
LL_bem = zeros(shz.N,src.N,3,3);
for component = 1:3
    for p = 1:numel(Kb)
        LL_bem(:,:,p,component) = LL_src_shz(:,:,p,component) + Kb{p}*bc.slip(:,:,component);
    end
end
end
