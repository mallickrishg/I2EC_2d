function matrixout = reconditionODEoperator(matrixinput)
% RECONDITIONODEOPERATOR  Stabilize an ODE operator matrix by removing spurious positive eigenvalues
%
%   matrixout = RECONDITIONODEOPERATOR(matrixinput)
%
%   This function takes a square matrix representing an ODE operator (e.g., a
%   stress kernel or decay matrix) and ensures numerical stability by:
%     1. Computing its eigenvalue decomposition.
%     2. Replacing any eigenvalues with non-negative real parts (>=0) by -eps.
%        This removes spurious growth modes that could cause instability in
%        exponential operators or periodic solutions.
%     3. Reconstructing the matrix from the corrected eigenvalues.
%
%   Input:
%     matrixinput  - (n x n) square matrix to be reconditioned.
%
%   Output:
%     matrixout    - (n x n) stabilized matrix with all eigenvalues strictly
%                    negative (or complex with negative real parts).
%
%   Notes:
%     - Only the real parts of the eigenvalues are checked; complex parts are
%       preserved.
%     - This is useful when small numerical errors produce spurious positive
%       eigenvalues that make exp(A*T) ill-conditioned or lead to unstable
%       periodic solves.
%
%   Example:
%     K_stable = reconditionODEoperator(K);
%
%   See also eig, expm, pinv

[Evectors,Evals] = eig(matrixinput);
% remove eigen values that cause instabilities 
Evals_corrected = diag(Evals);
Evals_corrected(real(Evals_corrected) >= 0) = -0.1;%-min(real(Evals_corrected))*1e5;
% reconstruct stress kernel
matrixout = real(Evectors * diag(Evals_corrected) / Evectors);

end