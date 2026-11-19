function y0 = computePeriodicInitialConditionBalanced(A, c, b, T)
% COMPUTEPERIODICINITIALCONDITIONBALANCED  Steady periodic initial state.
%
%   y0 = computePeriodicInitialConditionBalanced(A, c, b, T)
%
%   Solves the periodic offset condition y(0) = y(T) + b for the linear
%   ODE dy/dt = A*y + c on 0 <= t <= T: the initial state of an
%   earthquake cycle of period T whose coseismic jump in the state is b,
%   found without spinning up. Equivalently
%       (I - Phi)*y0 = eta + b,   Phi = expm(A*T),
%       eta = integral_0^T expm(A*(T-s))*c ds.
%   In modal coordinates (A = V*diag(lam)/V, zc = V\c, zb = V\b) each mode
%   is solved on its own:
%       z0 = (g(T).*zc + zb) ./ (1 - exp(lam*T)),
%   with g(T) from egh_modal and 1 - exp(lam*T) evaluated as -lam.*g(T).
%
% INPUTS
%   A       [n x n] ODE operator matrix (1/s)
%   c       [n x 1] constant forcing vector
%   b       [n x 1] periodic offset: coseismic jump in the state (1/s)
%   T       scalar period, positive (s)
%
% OUTPUTS
%   y0      [n x 1] initial state that satisfies y(0) = y(T) + b
%
% NOTES
%   - Stability: eigenvalues with Re(lam) > 0 are negated by
%     negateGrowingModes, exactly as in analyticSolveIVP, on the
%     assumption that they are spurious products of the discretisation.
%     Then |exp(lam*T)| < 1 for every mode with Re(lam) ~= 0, so the cycle
%     contracts and y0 is the state that a spin-up converges to. Build c
%     with the negated operator (see negateGrowingModes).
%   - Modes with |lam| < 1e-16 1/s (integrator modes that neither grow nor
%     decay) get z0 = 0, the minimum-norm solution. The threshold is
%     absolute, so the result is not continuous across it.
%   - Evaluating 1 - exp(lam*T) as -lam.*g(T) avoids cancellation for
%     small |lam*T|, where g(T) comes from a series.
%   - Warns if cond(V) > 1e14 (ill-conditioned eigenvectors).
%   - Changed in 2026: modes with Re(lam) > 1e-10 1/s were previously
%     dropped (z0 = 0) and modes with 0 < Re(lam) <= 1e-10 1/s solved with
%     their growing exponential. With negation the periodic state is
%     consistent with analyticSolveIVP (an IVP started from y0 returns to
%     y0 - b after one period) and depends continuously on the
%     viscosities, except where a complex-conjugate pair crosses
%     Re(lam) = 0 (see analyticSolveIVP).
%
% SEE ALSO
%   negateGrowingModes, egh_modal, analyticSolveIVP,
%   computePeriodicInitialCondition
%
% AUTHORS
%   Rishav Mallick (EOS), 2026; modified by Sharadha Sathiakumar (EOS), 2026

n = size(A, 1);

% ---- Eigendecomposition with growing modes negated ----
[V, lambda] = negateGrowingModes(A);

if cond(V) > 1e14
    warning('System remains ill-conditioned even after balancing.');
end

% ---- Modal coordinates ----
zc = V \ c;
zb = V \ b;

% ---- Periodic condition: z0 = (g(T)*zc + zb) / (1 - exp(lam*T)) ----
%  1 - exp(lam*T) = -lam * g(T) follows from g(T) = (exp(lam*T) - 1)/lam;
%  g(T) is accurate for all lam*T (series below |lam*T| = 1e-2, expm1
%  above), so the denominator is too.
[~, gT] = egh_modal(lambda, T);

z0 = zeros(n, 1);
active = abs(lambda) >= 1e-16;
one_minus_eT = -lambda(active) .* gT(active);
z0(active) = (gT(active) .* zc(active) + zb(active)) ./ one_minus_eT;

y0 = real(V * z0);

end
