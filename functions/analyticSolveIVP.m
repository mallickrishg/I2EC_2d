function [Y, U] = analyticSolveIVP(A, c, y0, tvec)
% ANALYTICSOLVEIVP  Solve the linear ODE dy/dt = A*y + c semi-analytically.
%
%   [Y, U] = analyticSolveIVP(A, c, y0, tvec)
%
%   Evaluates the state y(t) and its time integral u(t) at each time in
%   tvec. The operator is balanced and eigendecomposed with its growing
%   modes negated (negateGrowingModes); each mode is then evolved exactly
%   with the modal integrals of egh_modal:
%       z(t) = e(t).*z0 + g(t).*zc,    uz(t) = g(t).*z0 + h(t).*zc
%   with z0 = V\y0 and zc = V\c. This is the time stepper of the
%   semi-analytical (linear viscoelastic) path in I2EC.
%
% INPUTS
%   A       [n x n] ODE operator matrix (1/s), e.g. diag(1./eta)*K
%   c       [n x 1] constant forcing vector (units of dy/dt)
%   y0      [n x 1] initial condition at t = 0
%   tvec    [1 x m] evaluation times (s)
%
% OUTPUTS
%   Y       [n x m] state y(tvec(j)) in column j
%   U       [n x m] integral of y from 0 to tvec(j) in column j
%
% NOTES
%   - Stability: every eigenvalue with Re(lam) > 0 is negated
%     (lam -> -lam), on the assumption that such modes are spurious
%     products of the spatial discretisation (they depend on the mesh).
%     This is a regularisation, not a statement about the continuum
%     physics. No mode is dropped or frozen.
%   - Build the forcing c with the negated operator, e.g.
%     [~, ~, A] = negateGrowingModes(A); c = edot_lock - A*y_pl;
%     otherwise the long-term state y_pl is not a steady state.
%   - Outputs depend continuously on the viscosities except where a
%     complex-conjugate pair of eigenvalues crosses Re(lam) = 0: negation
%     flips the sign of its imaginary part, a jump of 2|sin(Im(lam) t)|
%     in that mode's amplitude (about 2|Im(lam)| t when small). A
%     spurious eigenvalue close to zero stays nearly neutral after
%     negation and decays slowly.
%   - Changed in 2026: modes with |Re(lam)| <= 1e-10*max(1,max|lam|) 1/s
%     were previously frozen (evolved as lam = 0) and modes with larger
%     Re(lam) > 0 dropped. Against the matrix exponential, the freeze
%     gave errors in mantle strain rate of up to 1.6% at 10 yr and 95%
%     over a 500-yr relaxation on the example mesh and three regenerated
%     meshes (497-1316 triangles; mantle 1e17-1e19 Pa s, uniform and
%     layered). The exact modal integrals agree with the matrix
%     exponential to round-off.
%   - Y and U are the real parts of the modal sums.
%
% SEE ALSO
%   negateGrowingModes, egh_modal, computePeriodicInitialConditionBalanced
%
% AUTHORS
%   Rishav Mallick (EOS), 2025; modified by Sharadha Sathiakumar (EOS), 2026

n = size(A, 1);

% ---- Eigendecomposition with growing modes negated ----
[V, lambda] = negateGrowingModes(A);

% ---- Modal coordinates ----
z0 = V \ y0;
ct = V \ c;

% ---- Evolve all modes with exact formulas ----
m = numel(tvec);
Z  = zeros(n, m);
UZ = zeros(n, m);

for j = 1:m
    [ej, gj, hj] = egh_modal(lambda, tvec(j));
    Z(:,j)  = ej .* z0 + gj .* ct;
    UZ(:,j) = gj .* z0 + hj .* ct;
end

% ---- Back to physical space ----
Y = real(V * Z);
U = real(V * UZ);

end
