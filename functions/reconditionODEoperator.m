function matrixout = reconditionODEoperator(matrixinput)
% RECONDITIONODEOPERATOR  Replace non-decaying eigenvalues by -0.1 (legacy).
%
%   matrixout = reconditionODEoperator(matrixinput)
%
%   Eigendecomposes a square ODE operator, sets every eigenvalue with a
%   non-negative real part to -0.1, and rebuilds the matrix:
%       matrixout = real(V*diag(lam_corrected)/V)
%
% INPUTS
%   matrixinput  [n x n] ODE operator matrix (1/s)
%
% OUTPUTS
%   matrixout    [n x n] reconditioned operator (1/s)
%
% NOTES
%   - Stability: a legacy damping rule, used only as the fallback in
%     computePeriodicInitialCondition. Eigenvalues with Re(lam) >= 0 are
%     replaced by the real number -0.1 1/s (a fixed, fast decay), so
%     their imaginary parts are lost. The current solvers use
%     negateGrowingModes (lam -> -lam for Re(lam) > 0) instead.
%   - No balancing is applied before eig, and the result is the real part
%     of the reconstruction.
%
% EXAMPLE
%   K_stable = reconditionODEoperator(K);
%
% SEE ALSO
%   negateGrowingModes, computePeriodicInitialCondition, eig
%
% AUTHORS
%   Rishav Mallick, 2025

[Evectors,Evals] = eig(matrixinput);
% remove eigen values that cause instabilities 
Evals_corrected = diag(Evals);
Evals_corrected(real(Evals_corrected) >= 0) = -0.1;%-min(real(Evals_corrected))*1e5;
% reconstruct stress kernel
matrixout = real(Evectors * diag(Evals_corrected) / Evectors);

end