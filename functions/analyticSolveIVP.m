function [Y, U] = analyticSolveIVP(A, c, y0, tvec)
% ANALYTICSOLVEIVP  Solve the linear ODE dy/dt = A*y + c semi-analytically.
%
%   [Y, U] = analyticSolveIVP(A, c, y0, tvec)
%
%   Evaluates the state y(t) and its time integral u(t) at each time in
%   tvec. The operator is balanced and eigendecomposed (eigBalanced); each
%   mode is then evolved exactly
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
%   - Stability: the operator must have no growing modes (Re(lam) >
%     1e-10 max|lam|), or the function raises an error
%     ('I2EC:growingModes'). The default kernels of
%     computeAllStressKernelsBem (cell-averaged, certified) guarantee this
%     for any viscosity. For centroid kernels ('cellaverage', 1) negate
%     the growing modes explicitly first and build the forcing with the
%     negated operator:
%         [~, ~, A] = negateGrowingModes(A); c = edot_lock - A*y_pl;
%     (Until 2026-11 this function negated growing modes itself.)
%   - With a negated operator, outputs jump where a complex-conjugate
%     pair of eigenvalues crosses Re(lam) = 0 as the viscosities change
%     (see negateGrowingModes).
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
%   eigBalanced, egh_modal, computePeriodicInitialConditionBalanced,
%   negateGrowingModes, computeAllStressKernelsBem
%
% AUTHORS
%   Rishav Mallick (EOS), 2025; modified by Sharadha Sathiakumar (EOS), 2026

n = size(A, 1);

% ---- Balanced eigendecomposition; growing modes are an error ----
[V, lambda, nGrowing] = eigBalanced(A);
if nGrowing > 0
    error('I2EC:growingModes', ['analyticSolveIVP: the operator has %d growing modes (largest ' ...
        'Re(lambda) = %.3g 1/s). Use the default kernels of computeAllStressKernelsBem ' ...
        '(cell-averaged, certified), or negate them explicitly with ' ...
        '[~, ~, A] = negateGrowingModes(A) and build the forcing with that A.'], ...
        nGrowing, max(real(lambda)));
end

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
