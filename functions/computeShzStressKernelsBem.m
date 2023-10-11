function LL = computeShzStressKernelsBem(src,shz,boundary)

% compute kernels relating unit shear from src to boundary
% TODO - create function to compute displacement kernel in half-space
LL_src_boundary = geometry.computeShzStressKernels(src,boundary);
G_src_boundary = geometry.computeShzDisplacementKernels(src,boundary.xc);

% compute displacement & traction kernels for boundary on itself
[Kdd,Kdn,Knd,Knn] = geometry.computeFullTractionKernels(boundary,boundary);
[Gdx,Gdz,Gnx,Gnz] = geometry.computeDisplacementKernels(boundary,boundary.xc);
Ktraction = [Kdd,Knd;Kdn,Knn];
Kdisp = [Gdx,Gnx;Gdz,Gnz];

% compute stress or traction kernels from shear zones to receiver 
% (could be shz or fault - check based on object type)

if isa(shz,'geometry.shearZoneReceiver')
    % TODO - create function to compute traction kernels from shz sources in a half-space
    LL_src_shz = geometry.computeShzStressKernels(src,shz);

    [Kdxx_boundary_shz,Kdzz_boundary_shz,Kdxz_boundary_shz,Knxx_boundary_shz,Knzz_boundary_shz,Knxz_boundary_shz] = ...
        geometry.computeFullStressKernels(boundary,shz);

elseif isa(shz,'geometry.receiver')
    
    LL_src_shz = geometry.computeShzStressKernels(src,shz);
    [Kdd_boundary_shz,Kdn_boundary_shz,Knd_boundary_shz,Knn_boundary_shz] = ...
        geometry.computeFullTractionKernels(boundary,shz);
else
    error('not a recognized geometry, provide either geometry.receiver or geometry.shearZoneReceiver')
end


end