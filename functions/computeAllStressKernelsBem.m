function evl = computeAllStressKernelsBem(rcv,shz,boundary,varargin)
% COMPUTEALLSTRESSKERNELSBEM  Assemble all stress interaction kernels.
%
%   evl = computeAllStressKernelsBem(rcv, shz, boundary)
%   evl = computeAllStressKernelsBem(..., 'cellaverage', 7)
%   evl = computeAllStressKernelsBem(..., 'cellaverage', 1, 'kernelcertify', false, 'kernelmodify', 1)
%
%   Default: mantle stress averaged over each cell with a 3-point rule and
%   the mantle-mantle kernel certified free of growing modes. The last
%   form gives the earlier centroid kernel with zeroed growing modes.
%
%   Computes the four blocks of the elastic interaction operator between
%   fault patches (rcv) and mantle shear-zone cells (shz), corrected for
%   the model boundary by the boundary element method: fault slip to fault
%   shear traction (KK), fault slip to mantle deviatoric stress (KL),
%   mantle deviatoric strain to fault shear traction (LK) and mantle
%   deviatoric strain to mantle deviatoric stress (LL). Mantle strain is
%   deviatoric (e33 = -e22), so each cell has two strain components
%   (e22, e23) and two stress components, sigma22' = (sxx - szz)/2 and
%   sigma23 = sxz.
%
% INPUTS
%   rcv       geometry.receiver object: megathrust patches (N patches)
%   shz       geometry.shearZoneReceiver object: mantle cells (M cells)
%   boundary  geometry.receiver object: model boundary patches
%             .Vpl  0: zero displacement; 1: zero traction (.Vx, .Vz
%             are ignored: the kernels solve the homogeneous problem)
%
% OPTIONS
%   cellaverage   stress in each mantle cell (KL and LL blocks) as
%                 1: the value at the centroid (collocation), or
%                 3, 7: the average over the cell with a 3- or 7-point
%                    interior quadrature rule  [3]
%   kernelcertify true: symmetrise the area-weighted LL and check that it
%                 is negative semidefinite (certifyMantleKernel), which
%                 guarantees no growing modes for any viscosity and
%                 power-law exponent; errors if the check fails. Needs
%                 'cellaverage' 3 or 7  [true]
%   kernelmodify  empirical fix for the growing modes of centroid
%                 kernels (see modifyGrowingModes); not needed with
%                 certified kernels  [0]
%                 0: keep LL as computed
%                 1: set the growing eigenvalues to 0 (neutral modes;
%                    they can accumulate over many cycles)
%                 2: set them to -max|lambda| (damped modes)
%                 1 and 2 are for numerical (ode45) solves only
%
% OUTPUTS
%   evl       struct with fields
%             .KK  [N x N] shear traction on rcv per unit slip on rcv
%                  (MPa/m)
%             .KL  [M x N x 2] deviatoric stress in shz per unit slip on
%                  rcv; (:,:,1) = (sxx - szz)/2, (:,:,2) = sxz (MPa/m)
%             .LK  [N x M x 2] shear traction on rcv per unit deviatoric
%                  strain in shz; (:,:,1) from e22 (with e33 = -e22),
%                  (:,:,2) from e23 (MPa)
%             .LL  [M x M x 2 x 2] deviatoric stress in shz per unit
%                  deviatoric strain in shz; 3rd index = receiver stress
%                  component (sigma22', sigma23), 4th index = source strain
%                  component (e22, e23) (MPa)
%             .certificate  result of certifyMantleKernel (certified
%                  kernels only)
%             .options  the cellaverage, kernelcertify and kernelmodify
%                  values used
%
% NOTES
%   - Cell averaging (default, 'cellaverage' 3) makes the mantle-mantle
%     block a Galerkin discretisation of the elastic operator (piecewise-
%     constant strain, cell-averaged stress), which keeps its symmetry and
%     sign. Centroid collocation ('cellaverage', 1) breaks both for
%     neighbouring cells and produces growing modes (43 on the example
%     mesh). Averaging costs about 3x (or 7x) the time of the KL and LL
%     blocks. KK and LK (fault receivers) are not averaged.
%   - Certification (default): the area-weighted LL is symmetrised and
%     checked to be negative semidefinite (certifyMantleKernel), which
%     guarantees no growing modes for any viscosity distribution and
%     power-law exponent; the call errors if the check fails. The result
%     is stored in evl.certificate.
%   - Kernel modification ('kernelmodify' 1 or 2, centroid kernels): the
%     2M x 2M matrix [LL11 LL12; LL21 LL22] is eigendecomposed, every
%     eigenvalue with positive real part is changed (modifyGrowingModes),
%     and the matrix is rebuilt (real part kept). Zeroed modes are neutral
%     and hold cell-scale noise in a spun-up cycle (strongly so when the
%     Maxwell time is much shorter than the recurrence); damped ones decay.
%     Kept as the earlier, empirical fix.
%   - KK, KL and LK are never modified.
%   - A comment in the code says the index construction of the deviatoric
%     LL kernel still needs to be checked.
%   - Units assume the shear modulus of the earth model is in MPa.
%   - Prints progress and timings to the command window.
%
% SEE ALSO
%   computeFaultTractionKernelsBem, computeShzStressKernelsBem,
%   computeAllDisplacementKernelsBem, assembleStressKernel,
%   modifyGrowingModes, negateGrowingModes
%
% AUTHORS
%   Rishav Mallick (Caltech) & Sharadha Sathiakumar (EOS), 2023

% Check the number of optional arguments
numOptionalArgs = nargin - 3; % Subtract the three required arguments

% default is to modify kernels
modifykernel = 0;
cellaverage = 3;
kernelcertify = true;
% Parse name-value pairs
for i = 1:2:numOptionalArgs
    argname = varargin{i};
    argvalue = varargin{i + 1};

    % Check the name and assign the value accordingly
    switch lower(argname)
        case 'kernelmodify'
            modifykernel = argvalue;
        case 'cellaverage'
            cellaverage = argvalue;
        case 'kernelcertify'
            kernelcertify = argvalue;        
        otherwise
            error(['Unknown parameter: ', argname]);
    end
end
if ~(isscalar(modifykernel) && ismember(modifykernel, [0 1 2]))
    error('kernelmodify must be 0, 1 or 2.');
end
if ~(isscalar(cellaverage) && ismember(cellaverage, [1 3 7]))
    error('cellaverage must be 1, 3 or 7.');
end
if ~(isscalar(kernelcertify) && (islogical(kernelcertify) || isnumeric(kernelcertify)))
    error('kernelcertify must be true or false.');
end
if kernelcertify && cellaverage == 1
    error(['kernelcertify needs cell-averaged kernels (''cellaverage'', 3 or 7): centroid kernels ' ...
           'are not negative semidefinite. For centroid kernels pass ''kernelcertify'', false.']);
end
[qpoints, qweights] = cellQuadrature(shz, cellaverage);
% the receiver side of the BEM kernels reads shz.xc only: evaluate at each
% quadrature point in turn and restore the centroids afterwards
xc0 = shz.xc;
restoreCentroids = onCleanup(@() setCentroids(shz, xc0));

% create data structure to hold all kernels
evl = [];

%% compute fault-fault interaction kernels
disp('Computing fault - fault traction kernels')

[K,~] = computeFaultTractionKernelsBem(rcv,rcv,boundary);

% store shear traction kernel in 'evl'
evl.KK = K;

%% compute fault-shz deviatoric interaction kernel

disp('Computing fault - shear zone traction kernels')
% note the ordering of kernels here is [sxx,szz,sxz]
LL1 = 0; LL2 = 0; LL3 = 0;
for q = 1:numel(qweights)
    shz.xc = qpoints{q};
    [K1,K3,K2] = computeFaultTractionKernelsBem(rcv,shz,boundary);
    LL1 = LL1 + qweights(q)*K1;
    LL2 = LL2 + qweights(q)*K2;
    LL3 = LL3 + qweights(q)*K3;
end
shz.xc = xc0;

% initialize deviatoric K-L kernel [Nshz x Nf x 2]
evl.KL = zeros(shz.N,rcv.N,2);
evl.KL(:,:,1) = (LL1-LL3)./2;
evl.KL(:,:,2) = LL2;

%% compute shz-shz interaction kernels

disp('Computing shear zone - shear zone stress kernels')
tic
LL = 0;
for q = 1:numel(qweights)
    shz.xc = qpoints{q};
    LL = LL + qweights(q)*computeShzStressKernelsBem(shz,shz,boundary);
end
shz.xc = xc0;

% deviatoric stress kernels [2x2]
evl.LL = zeros(shz.N,shz.N,2,2);

% construct deviatoric stress kernel and remove positive eigen values
% here the first 2 indices are for eigen strain source
%           last 2 indices are for the stress component
% evl.LL2333 corresponds to e23 source resulting in a change in s33
% in the large 3x3 matrix this would be in position [3-row,2-column]
% L2222 = (evl.LL2222 - evl.LL3322 - evl.LL2233 + evl.LL3333)./2;
% L2322 = (evl.LL2322 - evl.LL2333)./2;
% L2223 = (evl.LL2223 - evl.LL3323);
% L2323 = evl.LL2323;

% need to check/verify this construction (the indices always confuse me)
L11 = (LL(:,:,1,1) - LL(:,:,1,3) - LL(:,:,3,1) + LL(:,:,3,3))./2;
L12 = (LL(:,:,1,2) - LL(:,:,3,2))./2;
L21 = (LL(:,:,2,1) - LL(:,:,2,3));
L22 = LL(:,:,2,2);

% store deviatoric kernels in 'evl'
evl.LL(:,:,1,1) = L11;
evl.LL(:,:,1,2) = L12;
evl.LL(:,:,2,1) = L21;
evl.LL(:,:,2,2) = L22;

% cell-averaged kernels: symmetrise and certify negative semidefinite
if kernelcertify
    evl = certifyMantleKernel(evl, shz);
    fprintf('Mantle kernel certified: asymmetry removed %.2e, max eigenvalue %.2e of the spectrum, %d zero-energy modes clipped\n', ...
        evl.certificate.asymmetry, evl.certificate.maxEigenvalue, evl.certificate.clipped);
end

% change the growing modes (numerical solves only)
evl = modifyGrowingModes(evl, modifykernel);

% record how the kernels were built (stored by saveKernels)
evl.options = struct('cellaverage', cellaverage, 'kernelcertify', double(logical(kernelcertify)), ...
    'kernelmodify', modifykernel);
toc

%% compute stress interactions from shear zones to fault (LK kernels)

disp('Computing shear zone - fault traction kernels')
tic
LL = computeShzStressKernelsBem(shz,rcv,boundary);

% due to deviatoric state, e22 = -e33. Incorporating this into the kernels
% shape of kernel: [Nf x Nshz x 2]
evl.LK = zeros(rcv.N,shz.N,2);
% evl.LK(:,:,1) = LL(:,:,1) - LL(:,:,3);
% evl.LK(:,:,2) = LL(:,:,2);

evl.LK(:,:,1) = LL(:,:,1,1) - LL(:,:,1,3);
evl.LK(:,:,2) = LL(:,:,1,2);
toc

end

function [points, weights] = cellQuadrature(shz, order)
% interior quadrature points of every mantle cell, one [M x 2] array per
% point, and weights summing to 1
switch order
    case 1
        bary = [1 1 1]/3;
        weights = 1;
    case 3      % degree 2
        bary = [2/3 1/6 1/6; 1/6 2/3 1/6; 1/6 1/6 2/3];
        weights = [1 1 1]/3;
    case 7      % Dunavant degree 5
        a1 = 0.059715871789770; b1 = 0.470142064105115;
        a2 = 0.797426985353087; b2 = 0.101286507323456;
        bary = [1/3 1/3 1/3; a1 b1 b1; b1 a1 b1; b1 b1 a1; a2 b2 b2; b2 a2 b2; b2 b2 a2];
        weights = [0.225, 0.132394152788506*[1 1 1], 0.125939180544827*[1 1 1]];
end
points = cell(1, numel(weights));
for q = 1:numel(weights)
    points{q} = bary(q,1)*shz.A + bary(q,2)*shz.B + bary(q,3)*shz.C;
end
end

function setCentroids(shz, xc)
shz.xc = xc;
end
