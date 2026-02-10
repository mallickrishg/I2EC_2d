function y0 = computePeriodicInitialConditionBalanced(A, c, b, T)
%COMPUTEPERIODICINITIALCONDITION_BALANCED Compute periodic initial condition with automatic balancing
% 
% y0 = COMPUTEPERIODICINITIALCONDITION_BALANCED(A, c, b, T) solves for
% the initial condition of a linear ODE system with periodic offset constraint.
% 
% Solves the periodic offset condition:
%   y(0) = y(T) + b
% for the linear ODE:
%   dy/dt = A*y + c,    0 <= t <= T
% 
% The function computes:
%   Phi = expm(A*T)
%   eta = integral_0^T exp(A*(T-s)) * c ds
% and solves:
%   (I - Phi) * y0 = eta + b
% 
% This implementation uses automatic balancing and modal decomposition to
% handle stiff systems with widely varying timescales (e.g., fault zones
% with viscosity contrasts of 10^3 or more).
% 
% Inputs:
%   A - (n x n) ODE operator matrix (should have non-positive eigenvalues)
%   c - (n x 1) constant forcing vector
%   b - (n x 1) periodic offset vector
%   T - scalar, period duration (positive)
% 
% Output:
%   y0 - (n x 1) initial condition satisfying periodic constraint
% 
% Algorithm:
%   1. Apply diagonal balancing to A for better numerical conditioning
%   2. Compute eigendecomposition of balanced system
%   3. Solve in modal coordinates for each mode independently
%   4. Filter spurious unstable modes (real(lambda) > threshold)
%   5. Transform solution back to physical coordinates
% 
% Numerical Features:
%   - Handles nearly-singular (I - Phi) via modal decomposition
%   - Automatic filtering of numerically unstable eigenvalues
%   - Robust to viscosity contrasts exceeding 10^6
%   - Returns minimum-norm solution for near-zero eigenvalues
% 
% See also BALANCE, EIG, EXPM
% 
% Author: Rishav Mallick, EOS, 2026

n = size(A,1);

% Automated Balancing (Numerical Refinement)
[D_auto, Ab] = balance(A);

% Total transformation matrix: T_bal = S * D_auto
T_bal = D_auto;

% Eigendecomposition on the Balanced Matrix
[Vb, L] = eig(Ab);
lambda = diag(L);

% Recover Physical Eigenvectors
% If Ab = T_bal \ A * T_bal, then A's eigenvectors V = T_bal * Vb
V = T_bal * Vb;

% Check condition of corrected eigenvectors
if cond(V) > 1e14
    warning('System remains ill-conditioned even after balancing.');
end

% Proceed with your Modal Solver (using the refined V and lambda)
c_modal = V \ c;
b_modal = V \ b;

z0 = zeros(n, 1);
for i = 1:n
    % Apply your physics filters (Real(L) <= 0, etc.)
    if real(lambda(i)) <= 1e-10
        if abs(lambda(i)) < 1e-16
            z0(i) = 0; % Stable min-norm for integrators
        else
            phi_i = exp(lambda(i) * T);
            eta_i = ((phi_i - 1) / lambda(i)) * c_modal(i);
            z0(i) = (eta_i + b_modal(i)) / (1 - phi_i);
        end
    else
        z0(i) = 0; % Filter unstable physics
    end
end

y0 = real(V * z0);

end
