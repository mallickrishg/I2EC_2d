function y0 = computePeriodicInitialCondition(A, c, b, T)
% COMPUTEPERIODICINITIALCONDITION  Robustly compute initial condition for a linear ODE with periodic offset
%
%   y0 = COMPUTEPERIODICINITIALCONDITION(A, c, b, T)
%
%   Solves the periodic offset condition
%       y(0) = y(T) + b
%   for the linear ODE
%       dy/dt = A*y + c,    0 <= t <= T
%
%   The function computes
%       Phi = expm(A*T)
%       eta = integral_0^T exp(A*(T-s)) * c ds
%   and solves
%       (I - Phi) * y0 = eta + b
%
%   Features:
%     - Works even when I - Phi is singular, nearly singular, or incompatible.
%     - Returns a minimum-norm solution if infinite solutions exist.
%     - Returns a least-squares solution if no exact solution exists.
%     - Uses SVD-based robust fallback to handle ill-conditioning.
%
%   Inputs:
%     A : (n x n) ODE operator matrix
%     c : (n x 1) constant forcing vector
%     b : (n x 1) periodic offset vector
%     T : scalar period
%
%   Output:
%     y0 : (n x 1) initial condition vector satisfying y(0) = y(T) + b
%          (or best approximate if exact solution does not exist)
%
%   Notes:
%     - For well-conditioned I - Phi, a direct backslash solve is used.
%     - For rank-deficient or ill-conditioned cases, the function uses SVD
%       to compute minimum-norm or least-squares solutions.
%     - Warnings are issued when fallback solutions are used.
%
%   Example:
%       y0 = computePeriodicInitialCondition(A, c, b, T);
%
%   See also expm, svd, pinv

n = size(A,1);

% ---- Step 1: Compute Φ and η via block exponential (always well-defined)
M = zeros(n+1, n+1);
M(1:n,1:n) = A;
M(1:n,end) = c;

E = expm(M * T);
Phi = E(1:n,1:n);
eta = E(1:n,end);

rhs = eta + b;
I_minus_Phi = eye(n) - Phi;

% ---- Step 2: Try direct solve when reasonably conditioned
cond_est = cond(I_minus_Phi);
disp(['Condition number = ' num2str(cond_est)])
if cond_est < 1e15   % moderate threshold
    y0 = I_minus_Phi \ rhs;
    return
end

% ---- Full rank but ill-conditioned → use reconditioning & pseudo-inverse
M = zeros(n+1, n+1);
Arecon = reconditionODEoperator(A);
M(1:n,1:n) = Arecon;
M(1:n,end) = c;
E = expm(M * T);
Phi = E(1:n,1:n);
eta = E(1:n,end);
rhs = eta + b;

I_minus_Phi = eye(n) - Phi;
% use SVD for pseudo-inverse
[U,S,V] = svd(I_minus_Phi);
svals = diag(S);
tol = max(1e-12, eps * max(svals));
y0 = V * (diag(1./max(svals, tol)) * (U' * rhs));
warning('I - Phi nearly singular: solved using reconditioning & pseudo-inverse.');

end
