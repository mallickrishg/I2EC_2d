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
% 
% AUTHOR:
% Rishav Mallick, EOS, 2025

[V,D] = eig(A);
lambda = diag(D);

lam_scale = max(1,max(abs(lambda)));
eps_lam = 1e-10 * lam_scale;

stable   = real(lambda) < -eps_lam;
nearzero = abs(real(lambda)) <= eps_lam;
unstable = real(lambda) > eps_lam;

if any(unstable)
    warning('Discarding %d spurious unstable modes',sum(unstable))
end

% Modal coordinates
z0 = V \ y0;
ct = V \ c;

n = length(lambda);
m = numel(tvec);

Z = zeros(n,m);
UZ = zeros(n,m);

for j = 1:m
    t = tvec(j);

    % ---- stable modes ----
    ls = lambda(stable);
    es = exp(ls*t);

    Z(stable,j) = es .* z0(stable) ...
        + ((es-1)./ls) .* ct(stable);

    UZ(stable,j) = ((es-1)./ls) .* z0(stable) ...
        + ((es-1)./ls.^2 - t./ls) .* ct(stable);

    % ---- near-zero modes ----
    if any(nearzero)
        Z(nearzero,j)  = z0(nearzero) + t*ct(nearzero);
        UZ(nearzero,j) = t*z0(nearzero) + 0.5*t^2*ct(nearzero);
    end

    % ---- unstable modes: zero contribution ----
end

% Back to physical space
Y = real(V * Z);
U = real(V * UZ);

end