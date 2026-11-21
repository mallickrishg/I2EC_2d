function [V, lambda, nGrowing] = eigBalanced(A, tol, doBalance)
% EIGBALANCED  Eigendecomposition of an ODE operator (unbalanced by default).
%
%   [V, lambda] = eigBalanced(A)
%   [V, lambda, nGrowing] = eigBalanced(A, tol)
%   [V, lambda, nGrowing] = eigBalanced(A, tol, true)    balance first
%
%   Eigendecomposes A and counts its growing modes. Used by
%   analyticSolveIVP, computePeriodicInitialConditionBalanced and
%   negateGrowingModes, so all three see the same modes.
%
% INPUTS
%   A          [n x n] real ODE operator matrix (1/s)
%   tol        relative tolerance for a growing mode  [1e-10]
%   doBalance  true: balance A first (MATLAB balance) and return the
%              eigenvectors of A (V = T*Vb)  [false]
%
% OUTPUTS
%   V         [n x n] eigenvectors of A
%   lambda    [n x 1] eigenvalues of A (1/s)
%   nGrowing  number of eigenvalues with Re(lam) > tol*max|lam|. Smaller
%             positive real parts are taken as round-off.
%
% NOTES
%   - No balancing by default (since 2026-11-21). The I2EC operator is
%     A = K./eta: the kernel matrix, nearly symmetric once weighted by cell
%     area and patch width, with each row divided by a viscosity. Its
%     eigenvectors are then conditioned by the geometry, not by the
%     viscosities (cond(V) ~ 7e5 on the example for any mantle Maxwell
%     time). Balancing rescales rows and columns against that physical
%     scaling: with a weak fault and a strong mantle (Maxwell time 1e6 T)
%     cond(V) grew to 2e10 and the periodic mantle state lost three
%     digits (1e-6 instead of 1e-9 against a direct solve). The cost is
%     at moderate Maxwell times, where balancing was more accurate
%     (example 2: 2e-9 instead of 4e-12); unbalanced stays at 1e-9 or
%     better from t_M/T = 1e-5 to 1e9. A symmetrising diagonal scaling
%     (cond(V) ~ 10) loses digits for strong mantles like balancing, so
%     it is the relative accuracy of the small eigenvalues, kept by the
%     graded unbalanced matrix, that matters there.
%
% SEE ALSO
%   negateGrowingModes, analyticSolveIVP,
%   computePeriodicInitialConditionBalanced, balance, eig

if nargin < 2 || isempty(tol)
    tol = 1e-10;
end
if nargin < 3
    doBalance = false;
end
if doBalance
    [T_bal, Ab] = balance(A);
    [Vb, Db] = eig(Ab);
    V = T_bal * Vb;
else
    [V, Db] = eig(A, 'nobalance');
end
lambda = diag(Db);
nGrowing = sum(real(lambda) > tol * max(abs(lambda)));
end
