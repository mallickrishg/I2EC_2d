function y0 = computePeriodicInitialCondition(A, c, b, T)
% COMPUTEPERIODICINITIALCONDITION  Periodic initial state via expm (legacy).
%
%   y0 = computePeriodicInitialCondition(A, c, b, T)
%
%   Solves the periodic offset condition y(0) = y(T) + b for the linear
%   ODE dy/dt = A*y + c on 0 <= t <= T. Phi = expm(A*T) and
%   eta = integral_0^T expm(A*(T-s))*c ds come from one matrix exponential
%   of the augmented matrix [A c; 0 0], and
%       (I - Phi)*y0 = eta + b
%   is then solved directly. Superseded by
%   computePeriodicInitialConditionBalanced, which the drivers use.
%
% INPUTS
%   A       [n x n] ODE operator matrix (1/s)
%   c       [n x 1] constant forcing vector
%   b       [n x 1] periodic offset vector
%   T       scalar period (s)
%
% OUTPUTS
%   y0      [n x 1] initial state with y(0) = y(T) + b, or an approximate
%           solution when I - Phi is ill-conditioned
%
% NOTES
%   - Prints the condition number of I - Phi. If it is below 1e15, y0 is
%     found by a direct (backslash) solve.
%   - Otherwise A is replaced by reconditionODEoperator(A), which sets
%     every eigenvalue with Re(lam) >= 0 to -0.1 1/s; Phi and eta are
%     recomputed and y0 comes from an SVD pseudo-inverse with the
%     singular values floored at max(1e-12, eps*max(s)). A warning is
%     issued.
%   - Stability: unlike negateGrowingModes, this fallback replaces
%     non-decaying modes by a fixed decay rate of 0.1 1/s and drops their
%     imaginary parts, so it is not consistent with analyticSolveIVP.
%
% EXAMPLE
%   y0 = computePeriodicInitialCondition(A, c, b, T);
%
% SEE ALSO
%   computePeriodicInitialConditionBalanced, reconditionODEoperator, expm,
%   svd
%
% AUTHORS
%   Rishav Mallick, 2025

n = size(A,1);

% ---- Step 1: Compute Φ and η via block exponential (always well-defined)
M = zeros(n+1, n+1);
M(1:n,1:n) = A;
M(1:n,end) = c;

E = expm(M * T);
Phi = E(1:n,1:n);
eta = E(1:n,end);

rhs = eta + b;
I_minus_Phi = eye(n) - Phi;

% ---- Step 2: Try direct solve when reasonably conditioned
cond_est = cond(I_minus_Phi);
disp(['Condition number = ' num2str(cond_est)])
if cond_est < 1e15   % moderate threshold
    y0 = I_minus_Phi \ rhs;
    return
end

% ---- Full rank but ill-conditioned → use reconditioning & pseudo-inverse
M = zeros(n+1, n+1);
Arecon = reconditionODEoperator(A);
M(1:n,1:n) = Arecon;
M(1:n,end) = c;
E = expm(M * T);
Phi = E(1:n,1:n);
eta = E(1:n,end);
rhs = eta + b;

I_minus_Phi = eye(n) - Phi;
% use SVD for pseudo-inverse
[U,S,V] = svd(I_minus_Phi);
svals = diag(S);
tol = max(1e-12, eps * max(svals));
y0 = V * (diag(1./max(svals, tol)) * (U' * rhs));
warning('I - Phi nearly singular: solved using reconditioning & pseudo-inverse.');

end
