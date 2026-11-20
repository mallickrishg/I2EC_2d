function LO = computeShzDisplacementKernels(shz,obs)
% COMPUTESHZDISPLACEMENTKERNELS  Half-space displacement from zone strain.
%
%   LO = geometry.computeShzDisplacementKernels(shz, obs)
%
%   Displacement at the observation points due to a unit eigenstrain in
%   each triangle of shz, in a homogeneous elastic half-space. Each
%   triangle is loaded with unit exx, exz and ezz in turn, using
%   computeDisplacementPlaneStrainTriangleShearZone (Barbot, 2018).
%
% INPUTS
%   shz     geometry.shearZoneReceiver object: source triangles
%           .N           [1 x 1] number of triangles
%           .A, .B, .C   [M x 2] triangle vertices (m, z negative)
%           .earthModel  geometry.LDhs object (uses .nu)
%   obs     [Nobs x 2] observation points (x, z) (m, z <= 0)
%
% OUTPUTS
%   LO      [Nobs x M x 2 x 3] displacement per unit strain (m). The 3rd
%           index is the component (1: x, 2: z, positive up); the 4th is
%           the source strain (1: exx, 2: exz, 3: ezz). Example:
%           LO(:,:,1,3) is the x displacement due to ezz.
%
% NOTES
%   - Depths are made positive for the call, and the vertical
%     displacement is flipped back so that z is positive up.
%   - Uses parfor (parallel if a pool is available).
%
% SEE ALSO
%   computeDisplacementPlaneStrainTriangleShearZone,
%   computeShzDisplacementKernelsBem, geometry.computeShzStressKernels
%
% REFERENCES
%   Barbot, S. (2018). Bull. Seismol. Soc. Am., doi:10.1785/0120180058
%
% AUTHORS
%   Rishav Mallick (JPL), 2023

Nobs = length(obs(:,1));
nu = shz.earthModel.nu;

% initialize displacement kernels
LO1 = zeros(Nobs,shz.N);
LO2 = zeros(Nobs,shz.N);

% 2 x 3 displacemet kernels: 3 source strains and 2 components of
% displacement 
LO = zeros(Nobs,shz.N,2,3);

% source strain 100;010;001
I = eye(3);

tic
% convert all depths to positive numbers
obs(:,2) = -obs(:,2);

A = shz.A;A(:,2) = -A(:,2);
B = shz.B;B(:,2) = -B(:,2);
C = shz.C;C(:,2) = -C(:,2);

for i = 1:3
    % each iteration of 'i' goes through each eigen strain source
    % i = 1 corresponds to e22 source
    % i = 2 corresponds to e23 source
    % i = 3 corresponds to e33 source
    parfor k = 1:shz.N
        [u2,u3]=computeDisplacementPlaneStrainTriangleShearZone( ...
            obs(:,1),obs(:,2),...
            A(k,:),B(k,:),C(k,:),...
            I(i,1),I(i,2),I(i,3),...
            nu);

        % need to work with 2-d matrices because MATLAB doesn't like 3-d or
        % 4-d matrices inside parfor
        LO1(:,k) = u2(:);
        LO2(:,k) = -u3(:);
      
    end
    LO(:,:,1,i) = LO1; % horizontal component of displacement 
    LO(:,:,2,i) = LO2; % vertical component of displacement
end

end