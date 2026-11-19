function [V, lambda, A_neg] = negateGrowingModes(A)
% NEGATEGROWINGMODES  Eigendecompose A and negate its growing modes.
%
%   [V, lambda] = negateGrowingModes(A)
%   [V, lambda, A_neg] = negateGrowingModes(A)
%
%   Balances A (MATLAB balance), computes the eigendecomposition of the
%   balanced matrix and replaces every eigenvalue with Re(lam) > 0 by
%   -lam. The eigenvectors are kept, so
%       A_neg = V*diag(lambda)/V
%   has the same modes as A, with the growing ones turned into decaying
%   ones. analyticSolveIVP and computePeriodicInitialConditionBalanced
%   both call this function, so they treat eigenvalues identically.
%
% INPUTS
%   A       [n x n] real ODE operator matrix (1/s)
%
% OUTPUTS
%   V       [n x n] eigenvectors of A (V = T*Vb, with T from balance)
%   lambda  [n x 1] eigenvalues of A, with Re(lam) > 0 negated (1/s)
%   A_neg   [n x n] real negated operator V*diag(lambda)/V; computed only
%           when requested
%
% NOTES
%   - Stability: this is the damping rule of the semi-analytical path.
%     It assumes every eigenvalue with Re(lam) > 0 is spurious, i.e. a
%     product of the spatial discretisation (such eigenvalues depend on
%     the mesh). Negating them is a regularisation, not a claim about the
%     true continuum physics.
%   - A_neg is real because the eigenvalues stay in complex-conjugate
%     pairs after negation. Passing A_neg back to the solvers changes it
%     only by round-off.
%   - Build every forcing term that contains the operator with A_neg,
%     e.g. the long-term loading edot_lt = -A_neg*y_pl, so that the
%     plate-rate state y_pl stays the steady state of
%     dy/dt = A_neg*(y - y_pl) + ... With -A*y_pl instead,
%     (A_neg - A)*y_pl acts as an extra constant forcing on the negated
%     modes.
%
% SEE ALSO
%   analyticSolveIVP, computePeriodicInitialConditionBalanced, balance, eig
%
% AUTHORS
%   Sharadha Sathiakumar (EOS), 2026

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
