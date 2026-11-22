function [V, lambda, A_neg] = negateGrowingModes(A)
% NEGATEGROWINGMODES  Eigendecompose A and negate its growing modes.
%
%   [V, lambda] = negateGrowingModes(A)
%   [V, lambda, A_neg] = negateGrowingModes(A)
%
%   Computes the eigendecomposition of A (eigOperator) and replaces
%   every eigenvalue with Re(lam) > 0 by -lam. The eigenvectors are kept, so
%       A_neg = V*diag(lambda)/V
%   has the same modes as A, with the growing ones turned into decaying
%   ones.
%
%   Opt-in. The default kernels of computeAllStressKernelsBem are
%   cell-averaged and certified, so the operator has no growing modes and
%   needs no negation; solveAnalyticIvp and
%   computePeriodicInitialCondition do not negate and raise an
%   error if the operator has growing modes. For centroid kernels
%   ('cellaverage', 1) negate explicitly and pass A_neg to the solvers:
%       [~, ~, A] = negateGrowingModes(A);   % then build c with this A
%
% INPUTS
%   A       [n x n] real ODE operator matrix (1/s)
%
% OUTPUTS
%   V       [n x n] eigenvectors of A (from eigOperator)
%   lambda  [n x 1] eigenvalues of A, with Re(lam) > 0 negated (1/s)
%   A_neg   [n x n] real negated operator V*diag(lambda)/V; computed only
%           when requested
%
% NOTES
%   - Stability: this was the damping rule of the semi-analytical path
%     before certified kernels. It assumes every eigenvalue with
%     Re(lam) > 0 is spurious, i.e. a
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
%   eigOperator, solveAnalyticIvp, computePeriodicInitialCondition,
%   computeAllStressKernelsBem, modifyGrowingModes
%
% AUTHORS
%   Sharadha Sathiakumar (EOS), 2026

% ---- Eigendecomposition (A = V*diag(lambda)/V) ----
[V, lambda] = eigOperator(A);

% ---- Negate growing eigenvalues ----
pos = real(lambda) > 0;
lambda(pos) = -lambda(pos);

if nargout > 2
    % change only the negated modes: A_neg = A - 2*V(:,pos)*diag(lam0(pos))*W(pos,:)
    % (lam0: eigenvalues before negation; lambda(pos) = -lam0(pos) here)
    % with W = inv(V); with no growing mode A_neg is A itself, not a
    % reconstruction V*diag(lam)/V (which carries the eigendecomposition's
    % round-off, ~1e-9 in the mantle rows of the example operator)
    if any(pos)
        Winv = V \ eye(size(A));
        A_neg = A + real(V(:,pos) * diag(lambda(pos)) * Winv(pos,:)) * 2;
    else
        A_neg = A;
    end
end

end
