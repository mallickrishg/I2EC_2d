function gps = computeSiteVelocitiesBem(obs,rcv,shz,boundary,V,e22dot,e23dot)

% surface velocity 
% gps = [];
% gps.rcv.vh = (devl.KO(:,:,1)*V')'./Vpl;
% gps.rcv.vz = (devl.KO(:,:,2)*V')'./Vpl;
% 
% gps.shz.vh=(devl.LO(:,:,1,1)*e22dot' ...
%            +devl.LO(:,:,1,2)*e23dot')./Vpl;
% 
% gps.shz.vz=(devl.LO(:,:,2,1)*e22dot' ...
%            +devl.LO(:,:,2,2)*e23dot')./Vpl;

Nobs = length(obs(:,1));
Nt = length(V(:,1));

gps = [];
gps.vx = zeros(Nt,Nobs);
gps.vz = zeros(Nt,Nobs);

% compute kernels relating unit shear slip from rcv to boundary
[Kdd_rcv_boundary,Kdn_rcv_boundary,~,~] = geometry.computeFaultTractionKernels(rcv,boundary);
[Gdx_rcv_boundary,Gdz_rcv_boundary,~,~] = geometry.computeFaultDisplacementKernels(rcv,boundary.xc);

% compute kernels relating unit shear from shz to boundary
LL_shz_boundary = geometry.computeShzStressKernels(shz,boundary);
G_shz_boundary = geometry.computeShzDisplacementKernels(shz,boundary.xc);

% compute displacement & traction kernels for boundary on itself
[Kdd,Kdn,Knd,Knn] = geometry.computeFaultTractionKernels(boundary,boundary);
[Gdx,Gdz,Gnx,Gnz] = geometry.computeFaultDisplacementKernels(boundary,boundary.xc);
Ktraction = [Kdd,Knd;Kdn,Knn];
Kdisp = [Gdx,Gnx;Gdz,Gnz];

% src to rcv kernels
% and boundary to rcv kernels
[Gdx_rcv_obs,Gdz_rcv_obs,~,~] = geometry.computeFaultDisplacementKernels(rcv,obs);
G_shz_obs = geometry.computeShzDisplacementKernels(shz,obs);

[Gdx_boundary_obs,Gdz_boundary_obs,Gnx_boundary_obs,Gnz_boundary_obs] = ...
    geometry.computeFaultDisplacementKernels(boundary,obs);


for i = 1:Nt
    % calculate source tractions
    td_source = Kdd_rcv_boundary*V(i,:)' + ...
                LL_shz_boundary(:,:,1,1)*e22dot(i,:)' + ...
                LL_shz_boundary(:,:,1,2)*e23dot(i,:)' + ...
               -LL_shz_boundary(:,:,1,3)*e22dot(i,:)';
    tn_source = Kdn_rcv_boundary*V(i,:)' + ...
                LL_shz_boundary(:,:,2,1)*e22dot(i,:)' + ...
                LL_shz_boundary(:,:,2,2)*e23dot(i,:)' + ...
               -LL_shz_boundary(:,:,2,3)*e22dot(i,:)';
    % calculate source displacements
    ux_source = Gdx_rcv_boundary*V(i,:)' + ...
                G_shz_boundary(:,:,1,1)*e22dot(i,:)' + ...
                G_shz_boundary(:,:,1,2)*e23dot(i,:)' + ...
               -G_shz_boundary(:,:,1,3)*e22dot(i,:)';
    uz_source = Gdz_rcv_boundary*V(i,:)' + ...
                G_shz_boundary(:,:,2,1)*e22dot(i,:)' + ...
                G_shz_boundary(:,:,2,2)*e23dot(i,:)' + ...
               -G_shz_boundary(:,:,2,3)*e22dot(i,:)';

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
    disp = [Gdx_rcv_obs;Gdz_rcv_obs]*V(i,:)' + ...
           [G_shz_obs(:,:,1,1)-G_shz_obs(:,:,1,3),G_shz_obs(:,:,1,2);...
            G_shz_obs(:,:,2,1)-G_shz_obs(:,:,2,3),G_shz_obs(:,:,2,2)]*[e22dot(i,:)';e23dot(i,:)'] + ...
           [Gdx_boundary_obs,Gnx_boundary_obs;Gdz_boundary_obs,Gnz_boundary_obs]*bem_sol;
    
    gps.vx(i,:) = disp(1:end/2);
    gps.vz(i,:) = disp(end/2+1:end);

end   

end