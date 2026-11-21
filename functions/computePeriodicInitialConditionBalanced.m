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
%   - Stability: as in analyticSolveIVP, the operator must have no
%     growing modes, or the function raises an error ('I2EC:growingModes').
%     The default (cell-averaged, certified) kernels guarantee this; for
%     centroid kernels negate explicitly first (negateGrowingModes) and
%     build c with the negated operator. Then |exp(lam*T)| < 1 for every
%     mode with Re(lam) ~= 0, so the cycle contracts and y0 is the state
%     that a spin-up converges to.
%   - Neutral modes are an error too ('I2EC:neutralModes'): modes that
%     relax by less than 1e-12 of themselves over a cycle, |Re(lam)|*T <=
%     1e-12. Such a mode neither decays nor grows (an integrator, lam = 0,
%     or an undamped oscillation), so the periodic state is not
%     determined: it exists only if the forcing has no component along
%     the mode, and is then not unique. A spin-up would need more than
%     1e12 cycles. The test is on lam*T, not on lam relative to
%     max|lam|: a strong mantle with a weak fault has modes 1e-14 of
%     max|lam| that are still slow, real modes (Maxwell time 1e6 T:
%     min |lam|*T = 3e-9), and the modal solution is accurate for them
%     (agrees with a direct solve to 1e-6 in the mantle). The test is
%     also independent of the time unit.
%   - Evaluating 1 - exp(lam*T) as -lam.*g(T) avoids cancellation for
%     small |lam*T|, where g(T) comes from a series.
%   - Warns if cond(V) > 1e14 (ill-conditioned eigenvectors).
%   - The periodic state is consistent with analyticSolveIVP (an IVP
%     started from y0 returns to y0 - b after one period).
%
% SEE ALSO
%   eigBalanced, egh_modal, analyticSolveIVP, negateGrowingModes,
%   computePeriodicInitialCondition
%
% AUTHORS
%   Rishav Mallick (EOS), 2026; modified by Sharadha Sathiakumar (EOS), 2026


% ---- Balanced eigendecomposition; growing modes are an error ----
[V, lambda, nGrowing] = eigBalanced(A);
if nGrowing > 0
    error('I2EC:growingModes', ['computePeriodicInitialConditionBalanced: the operator has %d growing modes (largest ' ...
        'Re(lambda) = %.3g 1/s). Use the default kernels of computeAllStressKernelsBem ' ...
        '(cell-averaged, certified), or negate them explicitly with ' ...
        '[~, ~, A] = negateGrowingModes(A) and build the forcing with that A.'], ...
        nGrowing, max(real(lambda)));
end

tolNeutral = 1e-12;
neutral = abs(real(lambda)) * T <= tolNeutral;
if any(neutral)
    error('I2EC:neutralModes', ['computePeriodicInitialConditionBalanced: the operator has %d modes that ' ...
        'neither decay nor grow over a cycle (|Re(lambda)|*T <= %.0e, smallest |lambda| = %.3g 1/s), ' ...
        'so the periodic state is not determined. Check the kernel (computeAllStressKernelsBem ' ...
        'certifies it) and that every degree of freedom is restored by the operator.'], ...
        nnz(neutral), tolNeutral, min(abs(lambda)));
end

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

one_minus_eT = -lambda .* gT;
z0 = (gT .* zc + zb) ./ one_minus_eT;

y0 = real(V * z0);

end
