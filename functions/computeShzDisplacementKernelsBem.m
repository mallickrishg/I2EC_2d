function LL_bem = computeShzDisplacementKernelsBem(src,obs,boundary,scalar_value)
% COMPUTESHZDISPLACEMENTKERNELSBEM  Displacement kernels due to shear-zone strain, BEM-corrected.
%
%   LL_bem = computeShzDisplacementKernelsBem(src, obs, boundary, scalar_value)
%
%   For an eigenstrain of size scalar_value in each source triangle, one
%   strain component at a time, computes the half-space displacement at
%   the observation points and adds the response to the model boundary,
%   found by a boundary element solve.
%
% INPUTS
%   src           geometry.shearZoneReceiver object: source triangles (Ns)
%   obs           [Nobs x 2] observation points [x, z] (m)
%   boundary      geometry.receiver object: model boundary patches
%                 .Vpl  0: zero displacement; 1: zero traction
%   scalar_value  strain in each source triangle; 1 gives kernels per unit
%                 strain
%
% OUTPUTS
%   LL_bem        [Nobs x Ns x 2 x 3] displacement (m); 3rd index =
%                 displacement component (ux, uz), 4th index = source
%                 strain component (e22, e23, e33)
%
% NOTES
%   - Boundary correction as in computeFaultDisplacementKernelsBem (zero
%     displacement on Dirichlet patches; boundary.Vx, .Vz are not used).
%   - computeAllDisplacementKernelsBem combines the components into
%     deviatoric kernels (e33 = -e22).
%
% SEE ALSO
%   computeAllDisplacementKernelsBem, computeShzStressKernelsBem
%
% AUTHORS
%   Rishav Mallick, 2023

% these are [2 x 3] stress kernels
Nobs = length(obs(:,1));
LL_bem = zeros(Nobs,src.N,2,3);

% compute kernels relating unit shear from src to boundary
LL_src_boundary = geometry.computeShzStressKernels(src,boundary);
G_src_boundary = geometry.computeShzDisplacementKernels(src,boundary.xc);

% compute displacement & traction kernels for boundary on itself
[Kdd,Kdn,Knd,Knn] = geometry.computeFaultTractionKernels(boundary,boundary);
[Gdx,Gdz,Gnx,Gnz] = geometry.computeFaultDisplacementKernels(boundary,boundary.xc);
Ktraction = [Kdd,Knd;Kdn,Knn];
Kdisp = [Gdx,Gnx;Gdz,Gnz];

% compute stress or traction kernels from shear zone source to receiver 
LL_src_shz = geometry.computeShzDisplacementKernels(src,obs);

% compute kernels from boundary to shz receiver object
[Gdx_boundary_shz,Gdz_boundary_shz,Gnx_boundary_shz,Gnz_boundary_shz] = ...
        geometry.computeFaultDisplacementKernels(boundary,obs);


for component = 1:3
    Gx = zeros(Nobs,src.N);
    Gz = zeros(Nobs,src.N);
    for i = 1:src.N
        % we only consider shear sources
        source = zeros(src.N,1);
        source(i) = scalar_value;

        % calculate source tractions
        td_source = LL_src_boundary(:,:,1,component)*source;
        tn_source = LL_src_boundary(:,:,2,component)*source;
        % calculate source displacements
        ux_source = G_src_boundary(:,:,1,component)*source;
        uz_source = G_src_boundary(:,:,2,component)*source;

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
        
        % compute displacements from source term + bem solution on boundary
        disp = ...
            [LL_src_shz(:,:,1,component);...
            LL_src_shz(:,:,2,component)] * source + ...
            [Gdx_boundary_shz,Gnx_boundary_shz;...
            Gdz_boundary_shz,Gnz_boundary_shz] * bem_sol;

        Gx(:,i) = disp(1:end/2);% x-component
        Gz(:,i) = disp(end/2+1:end);% z-component
        
    end

    % store each K (kernel) in a bigger matrix
    LL_bem(:,:,1,component) = Gx; 
    LL_bem(:,:,2,component) = Gz; 

end