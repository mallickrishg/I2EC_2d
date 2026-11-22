function y0 = computePeriodicInitialCondition(A, c, b, T)
% COMPUTEPERIODICINITIALCONDITION  Steady periodic initial state.
%
%   y0 = computePeriodicInitialCondition(A, c, b, T)
%
%   Solves the periodic offset condition y(0) = y(T) + b for the linear
%   ODE dy/dt = A*y + c on 0 <= t <= T: the initial state of an
%   earthquake cycle of period T whose coseismic jump in the state is b,
%   found without spinning up. Equivalently
%       (I - Phi)*y0 = eta + b,   Phi = expm(A*T),
%       eta = integral_0^T expm(A*(T-s))*c ds.
%   Since Phi - I = A*G with G = integral_0^T expm(A*s) ds (so eta = G*c),
%   it is solved directly as
%       y0 = -(A*G) \ (G*c + b),
%   with G from one matrix exponential of the augmented matrix
%   [A I; 0 0]*T. The operator is still eigendecomposed (eigOperator), but
%   only to check for growing and neutral modes.
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
%   - Stability: as in solveAnalyticIvp, the operator must have no
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
%     min |lam|*T = 3e-9), for which the periodic state is still well
%     defined. The test is also independent of the time unit.
%   - Accuracy (changed 2026-11-22; was modal). The direct solve is
%     accurate for weak and strong mantles alike. On the example 2 operator
%     with mantle Maxwell times t_M/T = 1e-5 to 1e9 it agrees with the
%     balanced modal solve to 1e-11 or better where that one is accurate
%     (t_M/T <= 10) and with the unbalanced modal solve to ~1e-9 where
%     that one is (strong mantles); each modal variant alone lost up to
%     three digits (unbalanced: 1e-9 at t_M/T = 0.1, example 2: 3e-9) or
%     more (balanced: 8e-3 at 1e9). Forming I - Phi instead, with G*c
%     from an (n+1) exponential, loses digits for strong mantles.
%   - Cost: an expm of size 2n, about twice the eigendecomposition
%     (0.7 s against 0.3 s for n = 1061).
%   - The linear solve is equilibrated (MATLAB equilibrate): without it,
%     strong mantles (t_M/T >= 1e6) gave a "close to singular" warning
%     from the fault/mantle scale difference alone.
%   - The periodic state is consistent with solveAnalyticIvp (an IVP
%     started from y0 returns to y0 - b after one period).
%
% SEE ALSO
%   eigOperator, solveAnalyticIvp, negateGrowingModes, expm
%
% AUTHORS
%   Rishav Mallick (EOS), 2026; modified by Sharadha Sathiakumar (EOS), 2026


% ---- Eigendecomposition; growing modes are an error ----
[V, lambda, nGrowing] = eigOperator(A);
if nGrowing > 0
    error('I2EC:growingModes', ['computePeriodicInitialCondition: the operator has %d growing modes (largest ' ...
        'Re(lambda) = %.3g 1/s). Use the default kernels of computeAllStressKernelsBem ' ...
        '(cell-averaged, certified), or negate them explicitly with ' ...
        '[~, ~, A] = negateGrowingModes(A) and build the forcing with that A.'], ...
        nGrowing, max(real(lambda)));
end

tolNeutral = 1e-12;
neutral = abs(real(lambda)) * T <= tolNeutral;
if any(neutral)
    error('I2EC:neutralModes', ['computePeriodicInitialCondition: the operator has %d modes that ' ...
        'neither decay nor grow over a cycle (|Re(lambda)|*T <= %.0e, smallest |lambda| = %.3g 1/s), ' ...
        'so the periodic state is not determined. Check the kernel (computeAllStressKernelsBem ' ...
        'certifies it) and that every degree of freedom is restored by the operator.'], ...
        nnz(neutral), tolNeutral, min(abs(lambda)));
end

% ---- Periodic condition, solved directly ----
%  y(T) = Phi*y0 + G*c with Phi = expm(A*T), G = int_0^T expm(A*s) ds, and
%  y0 = y(T) + b give (I - Phi)*y0 = G*c + b. Since Phi - I = A*G exactly,
%      y0 = -(A*G) \ (G*c + b),
%  which avoids forming I - Phi (cancellation for slow modes). G comes from
%  one matrix exponential of the augmented 2n x 2n matrix [A I; 0 0]*T.
%  A*G inherits the ~1e6 scale difference between fault and mantle rows
%  (rcond down to 4e-19 for strong mantles, with MATLAB's "close to
%  singular" warning), so it is equilibrated (row and column scaling)
%  first: rcond ~1e-4 in every regime, same solution to 1e-15.
n = size(A, 1);
X = expm([A, eye(n); zeros(n, 2*n)] * T);
G = X(1:n, n+1:end);
M = A*G;
[P, R, C] = equilibrate(M);
y0 = -C * ((R*P*M*C) \ (R*P*(G*c + b)));

end
