function [Y, U] = analyticSolveIVP(A, c, y0, tvec)
% Computes y(t) and u(t) = ∫_0^t y(τ) dτ
%
% Uses augmented matrix exponential for the ODE system:
%   dy/dt = A y + c
%   du/dt = y
%
% Inputs:
%   A   : n×n matrix
%   c   : n×1 vector
%   y0  : n×1 vector (initial condition at t=0)
%   tvec: vector of times
%
% Outputs:
%   Y(:,j) = y(tvec(j))
%   U(:,j) = ∫_0^{tvec(j)} y(τ) dτ

n = size(A,1);
tvec = tvec(:).';
m = numel(tvec);

Y = zeros(n, m);
U = zeros(n, m);

% ---- Correct augmented matrix ----
% M = [ A   0   c
%       I   0   0
%       0   0   0 ]
M = zeros(2*n+1);
M(1:n,1:n)     = A;
M(1:n,2*n+1)   = c(:);
M(n+1:2*n,1:n) = eye(n);

% Initial augmented vector: [y0; u0; 1]
aug0 = [y0(:); zeros(n,1); 1];

for j = 1:m
    t = tvec(j);

    if t == 0
        Y(:,j) = y0;
        U(:,j) = 0;
        continue
    end

    E = expm(M * t);
    z = E * aug0;

    Y(:,j) = z(1:n);
    U(:,j) = z(n+1:2*n);
end

end