function LL = computeShzStressKernels(src,shz,method)
% COMPUTESHZSTRESSKERNELS  Half-space stress from strain in triangular zones.
%
%   LL = geometry.computeShzStressKernels(src, shz)
%   LL = geometry.computeShzStressKernels(src, shz, 'finitedifference')
%
%   Stress (or fault traction) at the receivers due to a unit eigenstrain
%   in each triangle of src, in a homogeneous elastic half-space. Each
%   source triangle is loaded with unit exx, exz and ezz, and the stress
%   is evaluated at the receiver centres with
%   computeStressPlaneStrainTriangleShearZone (closed-form derivative of
%   the displacement of Barbot, 2018), or with the older
%   computeStressPlaneStrainTriangleShearZoneFiniteDifference. If the
%   receivers are fault patches, the stress is projected to shear and
%   normal traction.
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
%   method  'analytic' (default) or 'finitedifference'
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
%   - 'analytic' differentiates the closed-form displacement exactly
%     (about 1e-12 relative error) and evaluates the three source strains
%     in one call, about 11 times faster than 'finitedifference', which
%     differences the displacement (about 1e-6 relative error).
%   - Depths are made positive (z -> -z) for the source vertices and
%     receiver points before the call, as that function expects.
%   - The component order is [sxx, sxz, szz]; geometry.computeFaultStressKernels
%     uses [sxx, szz, sxz].
%   - Uses parfor (parallel if a pool is available).
%
% SEE ALSO
%   computeStressPlaneStrainTriangleShearZone,
%   computeStressPlaneStrainTriangleShearZoneFiniteDifference,
%   computeShzStressKernelsBem, geometry.computeShzDisplacementKernels
%
% REFERENCES
%   Barbot, S. (2018). Bull. Seismol. Soc. Am., doi:10.1785/0120180058
%
% AUTHORS
%   Rishav Mallick (JPL), 2023

if nargin < 3
    method = 'analytic';
end
analytic = strcmpi(method, 'analytic');
if ~analytic && ~strcmpi(method, 'finitedifference')
    error('method must be ''analytic'' or ''finitedifference''.');
end

mu = src.earthModel.G;
nu = src.earthModel.nu;

% convert all depths to positive numbers
xc = shz.xc; xc(:,2) = -xc(:,2);

A = src.A;A(:,2) = -A(:,2);
B = src.B;B(:,2) = -B(:,2);
C = src.C;C(:,2) = -C(:,2);

isPatch = isa(shz,'geometry.receiver');
if isPatch
    nv = shz.nv; dv = shz.dv;
else
    nv = []; dv = [];
end

% K1, K2, K3 [Nr x M x 3]: the three stress (or traction) components;
% the 3rd index is the source strain (1: exx, 2: exz, 3: ezz)
K1 = zeros(shz.N,src.N,3);
K2 = zeros(shz.N,src.N,3);
K3 = zeros(shz.N,src.N,3);
I = eye(3);

parfor k = 1:src.N
    if analytic
        % all three unit source strains at once [Nr x 3]
        [s22,s23,s33] = computeStressPlaneStrainTriangleShearZone( ...
            xc(:,1),xc(:,2),A(k,:),B(k,:),C(k,:),I(1,:),I(2,:),I(3,:),mu,nu);
    else
        s22 = zeros(size(xc,1),3); s23 = s22; s33 = s22;
        for i = 1:3
            [s22(:,i),s23(:,i),s33(:,i)] = computeStressPlaneStrainTriangleShearZoneFiniteDifference( ...
                xc(:,1),xc(:,2),A(k,:),B(k,:),C(k,:),I(i,1),I(i,2),I(i,3),mu,nu);
        end
    end
    if isPatch
        % traction vector for the fault plane orientation, rotated to
        % fault-shear and fault-normal directions
        t1 = s22.*nv(:,1) - s23.*nv(:,2);
        t2 = -s23.*nv(:,1) + s33.*nv(:,2);
        K1(:,k,:) = reshape(dv(:,1).*t1 + dv(:,2).*t2, [], 1, 3);   % tau_shear
        K2(:,k,:) = reshape(nv(:,1).*t1 + nv(:,2).*t2, [], 1, 3);   % tau_normal
    else
        K1(:,k,:) = reshape(s22, [], 1, 3);    % sxx
        K2(:,k,:) = reshape(s23, [], 1, 3);    % sxz
        K3(:,k,:) = reshape(s33, [], 1, 3);    % szz
    end
end

LL = zeros(shz.N,src.N,3,3);
LL(:,:,1,:) = reshape(K1, shz.N, src.N, 1, 3);
LL(:,:,2,:) = reshape(K2, shz.N, src.N, 1, 3);
LL(:,:,3,:) = reshape(K3, shz.N, src.N, 1, 3);

end
