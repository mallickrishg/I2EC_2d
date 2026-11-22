function [V, lambda, nGrowing] = eigOperator(A, tol, doBalance)
% EIGOPERATOR  Eigendecomposition of an ODE operator, with or without balancing.
%
%   [V, lambda] = eigOperator(A)
%   [V, lambda, nGrowing] = eigOperator(A, tol)
%   [V, lambda, nGrowing] = eigOperator(A, tol, true)    balance first
%
%   Eigendecomposes A and counts its growing modes. Used by
%   solveAnalyticIvp, computePeriodicInitialCondition and
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
%   - Neither variant is accurate for every I2EC operator. A = K./eta has
%     fault rows ~1e6 larger than mantle rows. Without balancing the
%     mantle part of the modes loses digits for weak and moderate mantles
%     (example 2: 3e-9 in the evolution instead of 1e-13); with balancing
%     the slow modes of a strong mantle lose their relative accuracy
%     (Maxwell time 1e6 T and more: 1e-7 to 1e-2 in the periodic state).
%     So the callers do not rely on either alone (2026-11-22):
%     computePeriodicInitialCondition solves the periodic state directly
%     and uses this function only to check the modes, and
%     solveAnalyticIvp computes both variants and keeps the one whose
%     propagation of the initial state agrees better with expm.
%     negateGrowingModes uses the default (unbalanced).
%   - Fixed diagonal scalings (by the long-term state, by block, or a
%     symmetrising one) were tried and each fails in one regime too.
%
% SEE ALSO
%   negateGrowingModes, solveAnalyticIvp,
%   computePeriodicInitialCondition, balance, eig

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
