function F = shiftedModalIntegrals(lam, tau, k)
% SHIFTEDMODALINTEGRALS  Modal integrals minus their value at lam = 0.
%
%   F = shiftedModalIntegrals(lam, tau, k)
%
%   For every mode (rows) and segment length (columns), with x = lam*tau:
%     k = 1:  exp(x) - 1             = x*phi1(x)
%     k = 2:  g - tau                = tau*x*phi2(x)
%     k = 3:  h - tau^2/2            = tau^2*x*phi3(x)
%   where g, h are those of computeModalIntegrals and
%   phi_m(x) = sum_{q>=0} x^q/(q+m)!. Their divided differences in lam are
%   the pairwise integrals D_k (computePairIntegralSums), and subtracting
%   the lam = 0 value keeps the differences accurate for slow modes.
%   Series for |x| < 0.5 (24 terms), closed forms above.
lam = lam(:); tau = tau(:).';
x = lam .* tau;
small = abs(x) < 0.5;
switch k
    case 1
        F = expm1(x);
        return
    case 2
        phi = (expm1(x) - x) ./ x.^2;
        m = 2;
    case 3
        phi = (expm1(x) - x - x.^2/2) ./ x.^3;
        m = 3;
end
xs = x(small);
ps = zeros(size(xs)); f = factorial(m); xq = ones(size(xs));
for q = 0:23
    ps = ps + xq / f;
    xq = xq .* xs; f = f * (q + m + 1);
end
phi(small) = ps;
F = (tau.^(k - 1)) .* x .* phi;
end
