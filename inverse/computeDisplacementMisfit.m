function [J, res] = computeDisplacementMisfit(obs, fwd, dU)
% COMPUTEDISPLACEMENTMISFIT  Weighted misfit of surface-displacement time series, and its adjoint source.
%
%   [J, res] = computeDisplacementMisfit(obs, fwd)
%   [J, res] = computeDisplacementMisfit(obs, fwd, dU)
%
%   Observation layer of the adjoint path. Predicts the displacement
%       u(:,j) = obs.H*fwd.U(:,j) + obs.fixed(:,j) + o
%   at the observation times of solveLinearSequence, with an optional
%   nuisance offset o per data row (station and component) eliminated by
%   weighted least squares (variable projection), and returns
%       J = 1/2 * sum(W .* r.^2),   r = u - d,
%   and the adjoint source S = H'*(W.*r) for computeAdjointGradient.
%   Because the offsets are optimal, dJ/dtheta = S-based gradient
%   exactly (envelope theorem), and the Gauss-Newton Jacobian is the
%   tangent projected the same way.
%
%   Missing data: NaN in obs.d (or zero weight) marks an entry with no
%   observation, e.g. a station that starts recording midway, stops, or
%   has gaps. Missing entries carry no misfit and no adjoint source; each
%   row's offset is fitted to the entries it has, and a row with none
%   gets offset 0 and contributes nothing.
%
% INPUTS
%   obs    struct with fields
%          .H        [nd x n] station operator (assembleStationKernels)
%          .d        [nd x nObs] data, one column per fwd.seq.tObs;
%                    NaN where a station has no observation
%          .W        weights (inverse variances): scalar or [nd x nObs]
%                    [1]
%          .fixed    [nd x nObs] displacement that does not depend on
%                    viscosity (long-term trend, hinge, locked-patch
%                    slip)  [0]
%          .offsets  true: a free offset per data row  [true]
%   fwd    struct from solveLinearSequence
%   dU     optional [n x nObs x r] tangents (computeTangentSequence)
%
% OUTPUTS
%   J      misfit
%   res    struct with fields
%          .r        [nd x nObs] residuals (with the optimal offsets);
%                    0 where data are missing
%          .mask     [nd x nObs] true where an observation is used
%          .u        [nd x nObs] predicted displacement (with offsets)
%          .offset   [nd x 1] optimal offsets (0 without)
%          .S        [n x nObs] adjoint source dJ/dU
%          .rw       [nd*nObs x 1] whitened residual sqrt(W).*r
%          .Jw       [nd*nObs x r] whitened, offset-projected Jacobian
%                    of the residual along each tangent (with dU only):
%                    a Gauss-Newton step solves Jw*dxi = -rw
%
% SEE ALSO
%   assembleStationKernels, solveLinearSequence, computeAdjointGradient,
%   computeTangentSequence
%
% AUTHORS
%   Rishav Mallick (EOS), 2026

H = obs.H;
nObs = size(fwd.U, 2);
nd = size(H, 1);
W = 1; if isfield(obs, 'W') && ~isempty(obs.W), W = obs.W; end
W = W .* ones(nd, nObs);
d = obs.d;
mask = ~isnan(d) & W > 0;                     % observed entries
W(~mask) = 0;
d(~mask) = 0;
den = sum(W, 2);                              % per row; 0 for a row with no data
den(den == 0) = Inf;                          % its offset and projection are 0
fixed = 0; if isfield(obs, 'fixed') && ~isempty(obs.fixed), fixed = obs.fixed; end
useOffsets = ~isfield(obs, 'offsets') || isempty(obs.offsets) || obs.offsets;

u = H * fwd.U + fixed;
r = u - d;
offset = zeros(nd, 1);
if useOffsets
    offset = -sum(W .* r, 2) ./ den;
    r = r + offset;
end
r(~mask) = 0;
J = 0.5 * sum(W .* r.^2, 'all');
res.r = r;
res.u = u + offset;
res.offset = offset;
res.mask = mask;
res.S = H' * (W .* r);
res.rw = reshape(sqrt(W) .* r, [], 1);
if nargin > 2 && ~isempty(dU)
    nr = size(dU, 3);
    res.Jw = zeros(nd*nObs, nr);
    for k = 1:nr
        du = H * dU(:, :, k);
        if useOffsets
            du = du - sum(W .* du, 2) ./ den;
        end
        res.Jw(:, k) = reshape(sqrt(W) .* du, [], 1);
    end
end
end
