function [V, lambda, A_neg] = negateGrowingModes(A)
%NEGATEGROWINGMODES  Eigendecomposition of A with growing modes negated.
%
%   [V, lambda] = negateGrowingModes(A)
%   [V, lambda, A_neg] = negateGrowingModes(A)
%
% Balances A (MATLAB balance), computes its eigendecomposition and replaces
% every eigenvalue with Re(lam) > 0 by -lam.  The eigenvectors are kept, so
% A_neg = V diag(lambda) V^-1 has the same modes as A, with the growing
% ones turned into decaying ones.  A_neg is real because the eigenvalues
% stay in complex-conjugate pairs after negation.
%
% This assumes that every eigenvalue with Re(lam) > 0 is spurious: such
% eigenvalues arise from the spatial discretisation (they depend on the
% mesh).  Negating them is a regularisation, not a claim about the true
% continuum physics.
%
% analyticSolveIVP and computePeriodicInitialConditionBalanced call this
% function, so both treat eigenvalues identically.  Build any forcing
% term that contains the operator with A_neg, e.g. the long-term loading
%   edot_lt = -A_neg * y_pl,
% so that the long-term (plate-rate) state y_pl stays the steady state of
% the evolved system dy/dt = A_neg (y - y_pl) + ...  With -A * y_pl
% instead, the difference (A_neg - A) * y_pl acts as an extra constant
% forcing on the negated modes, and y_pl is no longer steady.
% Passing A_neg back to the solvers changes it only by round-off.
%
% Inputs:
%   A : (n x n) real ODE operator matrix
%
% Outputs:
%   V      : (n x n) eigenvectors of A (V = T * Vb, T from balance)
%   lambda : (n x 1) eigenvalues, with Re(lam) > 0 negated
%   A_neg  : (n x n) real negated operator V diag(lambda) V^-1
%            (computed only when requested)
%
% See also analyticSolveIVP, computePeriodicInitialConditionBalanced
%
% Author: Sharadha Sathiakumar, EOS, 2026

% ---- Balancing for numerical conditioning ----
[T_bal, Ab] = balance(A);

% ---- Eigendecomposition of balanced matrix ----
[Vb, Db] = eig(Ab);
lambda = diag(Db);

% ---- Recover physical eigenvectors: A = T_bal * Ab / T_bal ----
V = T_bal * Vb;

% ---- Negate growing eigenvalues ----
pos = real(lambda) > 0;
lambda(pos) = -lambda(pos);

if nargout > 2
    A_neg = real(V * diag(lambda) / V);
end

end
