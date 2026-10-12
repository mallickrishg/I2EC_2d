function y0 = computePeriodicInitialConditionBalanced(A, c, b, T)
%COMPUTEPERIODICINITIALCONDITIONBALANCED  Periodic initial condition with negate rule.
%
%   y0 = computePeriodicInitialConditionBalanced(A, c, b, T)
%
% Solves the periodic offset condition
%   y(0) = y(T) + b
% for the linear ODE  dy/dt = A y + c,  0 <= t <= T.
%
% Equivalently (I - Phi) y0 = eta + b, with Phi = expm(A*T) and
% eta = integral_0^T expm(A*(T-s)) c ds.  In modal coordinates
% (A = V diag(lam) V^-1, zc = V\c, zb = V\b) each mode is solved separately:
%   z0 = (g(T) * zc + zb) / (1 - exp(lam * T)),
% with g(T) = (exp(lam*T) - 1)/lam from egh_modal.  From that definition,
% 1 - exp(lam*T) = -lam * g(T), and the denominator is computed this way:
% g(T) is accurate also for small |lam*T| (series), so the denominator does
% not suffer from the cancellation in 1 - exp(lam*T).
%
% Eigenvalues with Re(lam) > 0 are negated (lam -> -lam) by
% negateGrowingModes, as in analyticSolveIVP, which assumes they are
% spurious (they arise from the spatial discretisation and depend on the
% mesh).  Then |exp(lam*T)| < 1 for every mode with Re(lam) ~= 0, so the
% cycle contracts and y0 is the state that a spin-up converges to.  Build
% the forcing c with the negated operator (see negateGrowingModes).
%
% What changed (2026) and why: previously, modes with Re(lam) > 1e-10 s^-1
% were dropped (z0 = 0) and modes with 0 < Re(lam) <= 1e-10 s^-1 were
% solved with their growing exponential.  With the negate rule, the
% periodic state is consistent with analyticSolveIVP (an IVP started from
% y0 returns to y0 - b after one period) and depends continuously on the
% viscosities, except where a complex-conjugate pair of eigenvalues
% crosses Re(lam) = 0 (negation flips the sign of its imaginary part; see
% analyticSolveIVP for the size of the jump).
% Negation is a regularisation of discretisation-induced growing modes,
% not a claim about the true continuum physics.  A spurious eigenvalue
% very close to zero stays nearly neutral after negation.
%
% Modes with |lam| < 1e-16 s^-1 are set to z0 = 0 (minimum-norm solution
% for integrator modes that neither grow nor decay), as before.  The
% threshold is absolute: a mode just above it is solved normally, so the
% result is not continuous across it.
%
% Inputs:
%   A - (n x n) ODE operator matrix
%   c - (n x 1) constant forcing vector
%   b - (n x 1) periodic offset (coseismic jump in strain rate)
%   T - scalar, period duration (positive, seconds)
%
% Output:
%   y0 - (n x 1) initial condition satisfying the periodic constraint
%
% See also negateGrowingModes, egh_modal, analyticSolveIVP
%
% Author: Rishav Mallick, EOS, 2026
% Modified: Sharadha Sathiakumar, EOS, 2026 (negate + exact solver)

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
