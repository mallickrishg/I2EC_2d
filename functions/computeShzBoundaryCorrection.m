function bc = computeShzBoundaryCorrection(src,boundary)
% COMPUTESHZBOUNDARYCORRECTION  Boundary slip that cancels each unit shear-zone strain on the model boundary.
%
%   bc = computeShzBoundaryCorrection(src, boundary)
%
%   For a unit eigenstrain in each source triangle, one strain component at
%   a time, solves for the slip on the boundary patches that makes the
%   total displacement zero on the Dirichlet patches and the total traction
%   zero on the Neumann patches. The result depends only on the sources
%   and the boundary, not on any receiver, so it can be computed once and
%   passed to computeShzStressKernelsBem for every set of receivers.
%
% INPUTS
%   src       geometry.shearZoneReceiver object: source triangles (Ns)
%   boundary  geometry.receiver object: model boundary patches (Nb)
%             .Vpl  0: zero displacement; 1: zero traction (.Vx, .Vz
%             are ignored: homogeneous problem)
%
% OUTPUTS
%   bc        struct with fields
%             .slip  [2Nb x Ns x 3] boundary slip per unit source strain;
%                    rows (dip slip; normal slip) of the boundary patches,
%                    3rd index = source strain component (e22, e23, e33)
%             .Ns, .Nb  numbers of sources and boundary patches (checked
%                    by computeShzStressKernelsBem)
%
% NOTES
%   - One LU factorisation of the 2Nb x 2Nb boundary matrix, solved for
%     all 3Ns right-hand sides at once.
%   - The source triangles enter through their vertices (.A, .B, .C) only,
%     so moving src.xc (cell quadrature) does not change the result.
%
% SEE ALSO
%   computeShzStressKernelsBem, computeAllStressKernelsBem,
%   computeFaultTractionKernelsBem
%
% AUTHORS
%   Rishav Mallick, 2026 (from the per-source solve of
%   computeShzStressKernelsBem, Rishav Mallick 2023)

% source to boundary: traction (dip, normal) and displacement (x, z)
LL_src_boundary = geometry.computeShzStressKernels(src,boundary);
G_src_boundary = geometry.computeShzDisplacementKernels(src,boundary.xc);

% boundary on itself
[Kdd,Kdn,Knd,Knn] = geometry.computeFaultTractionKernels(boundary,boundary);
[Gdx,Gdz,Gnx,Gnz] = geometry.computeFaultDisplacementKernels(boundary,boundary.xc);
Ktraction = [Kdd,Knd;Kdn,Knn];
Kdisp = [Gdx,Gnx;Gdz,Gnz];

% boundary conditions (only 2 types)
% boundary.Vpl = 0: dirichlet BC (u = 0)
% boundary.Vpl = 1: Neumann BC (t = 0)
idirichlet = boundary.Vpl(:) == 0;
ineumann = boundary.Vpl(:) == 1;
rowsD = [idirichlet;idirichlet];
rowsN = [ineumann;ineumann];

kernel = zeros(2*boundary.N,2*boundary.N);
kernel(rowsD,:) = Kdisp(rowsD,:);
kernel(rowsN,:) = Ktraction(rowsN,:);

% right-hand sides: minus the source displacement / traction, every
% source and strain component as a column
BC = zeros(2*boundary.N,src.N,3);
for component = 1:3
    ux = G_src_boundary(:,:,1,component);
    uz = G_src_boundary(:,:,2,component);
    td = LL_src_boundary(:,:,1,component);
    tn = LL_src_boundary(:,:,2,component);
    BC(rowsD,:,component) = -[ux(idirichlet,:);uz(idirichlet,:)];
    BC(rowsN,:,component) = -[td(ineumann,:);tn(ineumann,:)];
end

slip = kernel\reshape(BC,2*boundary.N,[]);
bc = struct('slip',reshape(slip,2*boundary.N,src.N,3),'Ns',src.N,'Nb',boundary.N);
end
