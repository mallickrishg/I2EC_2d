function [V, lambda, nGrowing] = eigBalanced(A, tol)
% EIGBALANCED  Balanced eigendecomposition of an ODE operator.
%
%   [V, lambda] = eigBalanced(A)
%   [V, lambda, nGrowing] = eigBalanced(A, tol)
%
%   Balances A (MATLAB balance), eigendecomposes the balanced matrix and
%   returns the eigenvectors of A (V = T*Vb) and its eigenvalues. Used by
%   analyticSolveIVP, computePeriodicInitialConditionBalanced and
%   negateGrowingModes, so all three see the same modes.
%
% INPUTS
%   A         [n x n] real ODE operator matrix (1/s)
%   tol       relative tolerance for a growing mode  [1e-10]
%
% OUTPUTS
%   V         [n x n] eigenvectors of A
%   lambda    [n x 1] eigenvalues of A (1/s)
%   nGrowing  number of eigenvalues with Re(lam) > tol*max|lam|. Smaller
%             positive real parts are round-off of zero eigenvalues
%             (e.g. the zero-energy modes of a certified kernel).
%
% SEE ALSO
%   negateGrowingModes, analyticSolveIVP,
%   computePeriodicInitialConditionBalanced, balance, eig

if nargin < 2
    tol = 1e-10;
end
[T_bal, Ab] = balance(A);
[Vb, Db] = eig(Ab);
lambda = diag(Db);
V = T_bal * Vb;
nGrowing = sum(real(lambda) > tol * max(abs(lambda)));
end
