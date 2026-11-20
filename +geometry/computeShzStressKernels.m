function LL = computeShzStressKernels(src,shz)
% COMPUTESHZSTRESSKERNELS  Half-space stress from strain in triangular zones.
%
%   LL = geometry.computeShzStressKernels(src, shz)
%
%   Stress (or fault traction) at the receivers due to a unit eigenstrain
%   in each triangle of src, in a homogeneous elastic half-space. Each
%   source triangle is loaded with unit exx, exz and ezz in turn, and the
%   stress is evaluated at the receiver centres with
%   computeStressPlaneStrainTriangleShearZoneFiniteDifference (Barbot,
%   2018). If the receivers are fault patches, the stress is projected to
%   shear and normal traction.
%
% INPUTS
%   src     geometry.shearZoneReceiver object: source triangles
%           .N           [1 x 1] number of triangles
%           .A, .B, .C   [M x 2] triangle vertices (m, z negative)
%           .earthModel  geometry.LDhs object (uses .G, .nu)
%   shz     receivers: geometry.shearZoneReceiver (stress at the triangle
%           centroids) or geometry.receiver (traction on the patches)
%           .N           [1 x 1] number of receivers
%           .xc          [Nr x 2] receiver points (m, z negative)
%           .dv, .nv     [Nr x 2] unit dip and normal vectors (patches
%                        only)
%
% OUTPUTS
%   LL      [Nr x M x 3 x 3] kernel (MPa per unit strain). The 4th index is
%           the source strain (1: exx, 2: exz, 3: ezz). The 3rd index is
%           the stress component (1: sxx, 2: sxz, 3: szz) for
%           shearZoneReceiver receivers, or (1: shear traction, 2: normal
%           traction, 3: zero) for receiver patches. Example: LL(:,:,3,2)
%           is szz due to exz.
%
% NOTES
%   - Stress comes from finite-difference derivatives of the analytical
%     displacement field, not from a closed-form stress solution (see
%     computeStressPlaneStrainTriangleShearZoneFiniteDifference).
%   - Depths are made positive (z -> -z) for the source vertices and
%     receiver points before the call, as that function expects.
%   - The component order is [sxx, sxz, szz]; geometry.computeFaultStressKernels
%     uses [sxx, szz, sxz].
%   - Uses parfor (parallel if a pool is available).
%
% SEE ALSO
%   computeStressPlaneStrainTriangleShearZoneFiniteDifference,
%   computeShzStressKernelsBem, geometry.computeShzDisplacementKernels
%
% REFERENCES
%   Barbot, S. (2018). Bull. Seismol. Soc. Am., doi:10.1785/0120180058
%
% AUTHORS
%   Rishav Mallick (JPL), 2023

mu = src.earthModel.G;
nu = src.earthModel.nu;

% initialize stress kernels
LL1 = zeros(shz.N,src.N);
LL2 = zeros(shz.N,src.N);
LL3 = zeros(shz.N,src.N);

% this are 3x3 stress kernels
LL = zeros(shz.N,src.N,3,3);

% source strain 100;010;001
I = eye(3);

% convert all depths to positive numbers
xc = shz.xc; xc(:,2) = -xc(:,2);

A = src.A;A(:,2) = -A(:,2);
B = src.B;B(:,2) = -B(:,2);
C = src.C;C(:,2) = -C(:,2);

for i = 1:3
    % each iteration of 'i' goes through each eigen strain source
    % i = 1 corresponds to exx source
    % i = 2 corresponds to exz source
    % i = 3 corresponds to ezz source
    parfor k = 1:src.N
        [s22,s23,s33] = computeStressPlaneStrainTriangleShearZoneFiniteDifference( ...
            xc(:,1),xc(:,2),...
            A(k,:),B(k,:),C(k,:),...
            I(i,1),I(i,2),I(i,3),...
            mu,nu);
        
        % need to work with 2-d matrices because MATLAB doesn't like 3-d or
        % 4-d matrices inside parfor
        if isa(shz,'geometry.shearZoneReceiver')
            LL1(:,k) = s22(:);
            LL2(:,k) = s23(:);
            LL3(:,k) = s33(:);
        elseif isa(shz,'geometry.receiver')
            % compute traction vector for fault plane orientation
            t=[s22.*shz.nv(:,1) - s23.*shz.nv(:,2), ...
                -s23.*shz.nv(:,1) + s33.*shz.nv(:,2)];
            % rotate traction vector to fault-shear & fault-normal direction
            LL1(:,k) = shz.dv(:,1).*t(:,1) + shz.dv(:,2).*t(:,2);
            LL2(:,k) = shz.nv(:,1).*t(:,1) + shz.nv(:,2).*t(:,2);
        else
            LL1(:,k) = s22(:);
            LL2(:,k) = s23(:);
            LL3(:,k) = s33(:);
            % error('not a recognized geometry. provide either geometry.shearZoneReceiver or geometry.receiver')
        end

    end
    if isa(shz,'geometry.shearZoneReceiver')
        LL(:,:,1,i) = LL1; % sxx
        LL(:,:,2,i) = LL2; % sxz
        LL(:,:,3,i) = LL3; % szz
    elseif isa(shz,'geometry.receiver')
        LL(:,:,1,i) = LL1; % tau_shear
        LL(:,:,2,i) = LL2; % tau_normal
    else
        LL(:,:,1,i) = LL1; % sxx
        LL(:,:,2,i) = LL2; % sxz
        LL(:,:,3,i) = LL3; % szz
        % error('not a recognized geometry. provide either geometry.shearZoneReceiver or geometry.receiver')
    end
end

end