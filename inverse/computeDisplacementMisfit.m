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
% INPUTS
%   obs    struct with fields
%          .H        [nd x n] station operator (assembleStationKernels)
%          .d        [nd x nObs] data, one column per fwd.seq.tObs
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
%          .r        [nd x nObs] residuals (with the optimal offsets)
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
fixed = 0; if isfield(obs, 'fixed') && ~isempty(obs.fixed), fixed = obs.fixed; end
useOffsets = ~isfield(obs, 'offsets') || isempty(obs.offsets) || obs.offsets;

u = H * fwd.U + fixed;
r = u - obs.d;
offset = zeros(nd, 1);
if useOffsets
    offset = -sum(W .* r, 2) ./ sum(W, 2);
    r = r + offset;
end
J = 0.5 * sum(W .* r.^2, 'all');
res.r = r;
res.u = u + offset;
res.offset = offset;
res.S = H' * (W .* r);
res.rw = reshape(sqrt(W) .* r, [], 1);
if nargin > 2 && ~isempty(dU)
    nr = size(dU, 3);
    res.Jw = zeros(nd*nObs, nr);
    for k = 1:nr
        du = H * dU(:, :, k);
        if useOffsets
            du = du - sum(W .* du, 2) ./ sum(W, 2);
        end
        res.Jw(:, k) = reshape(sqrt(W) .* du, [], 1);
    end
end
end
