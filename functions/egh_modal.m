function [e, g, h] = egh_modal(lam, t)
%EGH_MODAL  Exact modal integrals for the linear ODE  dy/dt = A y + c.
%
%   [e, g, h] = egh_modal(lam, t)
%
% For each eigenvalue lam and scalar time t, returns
%   e = exp(lam * t)
%   g = integral_0^t exp(lam * s) ds        = (exp(lam*t) - 1) / lam
%   h = integral_0^t g(s) ds                = ((exp(lam*t) - 1)/lam - t) / lam
%
% All three are continuous through lam = 0 (g -> t, h -> t^2/2).  With
% x = lam * t, an 8-term Taylor series is used when |x| < 1e-2 (|x|
% includes the imaginary part; truncation error < 1e-21 at |x| = 1e-2),
% and the closed forms g = expm1(x)/lam, h = (expm1(x) - x)/lam^2 above.
% expm1 avoids the cancellation in exp(x) - 1; h still cancels for small
% |x| (relative error ~ eps/|x|), hence the series below the switch.
%
% Inputs:
%   lam : (n x 1) complex eigenvalues
%   t   : scalar, time (seconds)
%
% Outputs:
%   e, g, h : (n x 1) complex vectors
%
% Used by analyticSolveIVP and computePeriodicInitialConditionBalanced
% (eigenvalues from negateGrowingModes).
%
% Author: Sharadha Sathiakumar, EOS, 2026

x = lam * t;
e = exp(x);
g = zeros(size(lam));
h = zeros(size(lam));

s = abs(x) < 1e-2;

% ---- Series branch (|lam*t| < 1e-2): 8 terms ----
%   g/t   = sum_{k=0}^{7} x^k / (k+1)!
%   h/t^2 = sum_{k=0}^{7} x^k / ((k+1)! * (k+2))
xs = x(s);
gs = zeros(size(xs));
hs = zeros(size(xs));
fk = 1;                                % fk = (k+1)!
for k = 0:7
    gs = gs + xs.^k / fk;
    hs = hs + xs.^k / (fk * (k + 2));
    fk = fk * (k + 2);
end
g(s) = t * gs;
h(s) = t^2 * hs;

% ---- Closed-form branch (|lam*t| >= 1e-2) ----
xb = x(~s);
em1 = expm1(xb);
g(~s) = em1 ./ lam(~s);
h(~s) = (em1 - xb) ./ lam(~s).^2;

end
