function [Y, U] = analyticSolveIVP(A, c, y0, tvec)
%ANALYTICSOLVEIVP  Solve the linear ODE  dy/dt = A y + c  semi-analytically.
%
%   [Y, U] = analyticSolveIVP(A, c, y0, tvec)
%
% Computes y(t) and u(t) = integral_0^t y(s) ds at each time in tvec by
% eigendecomposition of the balanced operator and exact modal integrals
% (egh_modal).
%
% Eigenvalues with Re(lam) > 0 are negated (lam -> -lam) by
% negateGrowingModes, which assumes they are spurious (they arise from the
% spatial discretisation and depend on the mesh).  Every mode is then
% evolved exactly: no mode is dropped or frozen.  Build the forcing c with
% the negated operator (see negateGrowingModes), e.g.
%   [~, ~, A] = negateGrowingModes(A);   c = edot_lock - A * y_pl;
%
% What changed (2026) and why:
%   - Previously, modes with |Re(lam)| <= 1e-10*max(1,max|lam|) s^-1 were
%     frozen (evolved as lam = 0) and modes with larger Re(lam) > 0 were
%     dropped, with a warning.  The freeze was the solver's only error.
%     In a forced relaxation from the plate-rate state plus a coseismic
%     step, on the example mesh and three regenerated meshes (497-1316
%     triangles; uniform mantle
%     1e17-1e19 Pa s and layered), the error in mantle strain rate,
%     ||error|| / ||reference|| over all mantle cells against the matrix
%     exponential, was up to 1.6% at 10 yr after the earthquake and up to
%     95% over the 500-yr relaxation; on meshes with no mode in the freeze
%     band it was <= 5e-9.  The exact integrals remove it (round-off
%     agreement with the matrix exponential).
%   - Negating instead of dropping gives a contracting operator, a
%     consistent periodic cycle (computePeriodicInitialConditionBalanced
%     uses the same rule) and outputs that depend continuously on the
%     viscosities, except where a complex-conjugate pair of eigenvalues
%     crosses Re(lam) = 0: there negation flips the sign of its imaginary
%     part, so that mode's response changes from exp(i Im(lam) t) to
%     exp(-i Im(lam) t), a jump of 2|sin(Im(lam) t)| of its amplitude
%     (about 2|Im(lam)| t when that is small).  Negation is a
%     regularisation of discretisation-induced growing modes, not a claim
%     about the true continuum physics.  A spurious eigenvalue very close
%     to zero stays nearly neutral after negation, so it decays slowly.
%
% Inputs:
%   A    : (n x n) ODE operator matrix
%   c    : (n x 1) constant forcing vector
%   y0   : (n x 1) initial condition at t = 0
%   tvec : (1 x m) vector of evaluation times (seconds)
%
% Outputs:
%   Y(:,j) = y(tvec(j))                     state at each time
%   U(:,j) = integral_0^{tvec(j)} y(s) ds   integrated state
%
% See also negateGrowingModes, egh_modal, computePeriodicInitialConditionBalanced
%
% AUTHOR:
% Rishav Mallick, EOS, 2025
% Modified: Sharadha Sathiakumar, EOS, 2026 (negate + exact solver)

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
