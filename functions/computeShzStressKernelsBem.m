function LL_bem = computeShzStressKernelsBem(src,shz,boundary)

% these are [3 x 3] or [2 x 3] stress kernels
LL_bem = zeros(shz.N,src.N,3,3);

% compute kernels relating unit shear from src to boundary
LL_src_boundary = geometry.computeShzStressKernels(src,boundary);
G_src_boundary = geometry.computeShzDisplacementKernels(src,boundary.xc);

% compute displacement & traction kernels for boundary on itself
[Kdd,Kdn,Knd,Knn] = geometry.computeFaultTractionKernels(boundary,boundary);
[Gdx,Gdz,Gnx,Gnz] = geometry.computeFaultDisplacementKernels(boundary,boundary.xc);
Ktraction = [Kdd,Knd;Kdn,Knn];
Kdisp = [Gdx,Gnx;Gdz,Gnz];

% compute stress or traction kernels from shear zone source to receiver 
LL_src_shz = geometry.computeShzStressKernels(src,shz);

% compute kernels from boundary to shz receiver object
% (receiver could be shz or fault - check based on object type)
if isa(shz,'geometry.shearZoneReceiver')    
    [Kdxx_boundary_shz,Kdzz_boundary_shz,Kdxz_boundary_shz,Knxx_boundary_shz,Knzz_boundary_shz,Knxz_boundary_shz] = ...
        geometry.computeFaultStressKernels(boundary,shz);

elseif isa(shz,'geometry.receiver')
    [Kdd_boundary_shz,Kdn_boundary_shz,Knd_boundary_shz,Knn_boundary_shz] = ...
        geometry.computeFaultTractionKernels(boundary,shz);
else
    error('not a recognized geometry, provide either geometry.receiver or geometry.shearZoneReceiver')
end

for component = 1:3
    Kxx = zeros(shz.N,src.N);
    Kzz = zeros(shz.N,src.N);
    Kxz = zeros(shz.N,src.N);
    for i = 1:src.N
        % we only consider shear sources
        source = zeros(src.N,1);
        source(i) = 1;

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
        BC([idirichlet;idirichlet]) = [boundary.Vx(idirichlet);boundary.Vz(idirichlet)] - ...
                                      [ux_source(idirichlet);uz_source(idirichlet)];
        % traction BC
        BC([ineumann;ineumann]) = 0 - [td_source(ineumann);tn_source(ineumann)];

        % construct BEM kernel
        kernel = zeros(boundary.N*2,boundary.N*2);
        kernel([idirichlet;idirichlet],:) = Kdisp([idirichlet;idirichlet],:);
        kernel([ineumann;ineumann],:) = Ktraction([ineumann;ineumann],:);

        % BEM solve
        bem_sol = kernel\BC;
        
        % compute tractions or stress from source term + bem solution on boundary
        if isa(shz,'geometry.shearZoneReceiver')
            sigma = ...
                [LL_src_shz(:,:,1,component);...
                LL_src_shz(:,:,2,component);...
                LL_src_shz(:,:,3,component)] * source + ...
                [Kdxx_boundary_shz,Knxx_boundary_shz;...
                Kdzz_boundary_shz,Knzz_boundary_shz;...
                Kdxz_boundary_shz,Knxz_boundary_shz] * bem_sol;

            Kxx(:,i) = sigma(1:end/3);
            Kxz(:,i) = sigma(end/3+1:2/3*end);
            Kzz(:,i) = sigma(2/3*end+1:end);
        elseif isa(shz,'geometry.receiver')
            tau = ...
                [LL_src_shz(:,:,1,component);...
                LL_src_shz(:,:,2,component)] * source + ...
                [Kdd_boundary_shz,Knd_boundary_shz;...
                Kdn_boundary_shz,Knn_boundary_shz] * bem_sol;

            Kxz(:,i) = tau(1:end/2);% tau_shear
            Kxx(:,i) = tau(end/2+1:end);% tau_normal
        else
            error('not a recognized geometry, provide either geometry.receiver or geometry.shearZoneReceiver')
        end
    end

    % store each K (kernel) in a bigger matrix
    if isa(shz,'geometry.shearZoneReceiver')
        LL_bem(:,:,1,component) = Kxx; % sxx
        LL_bem(:,:,2,component) = Kxz; % sxz
        LL_bem(:,:,3,component) = Kzz; % szz
    elseif isa(shz,'geometry.receiver')
        LL_bem(:,:,1,component) = Kxz; % tau_shear
        LL_bem(:,:,2,component) = Kxx; % tau_normal
    else
        error('not a recognized geometry. provide either geometry.shearZoneReceiver or geometry.receiver')
    end
end

