function [D1, D2, D3] = computeExpDividedDifferences(a, b, tau)
% COMPUTEEXPDIVIDEDDIFFERENCES  Divided differences of exp(x*tau) on node pairs.
%
%   [D1, D2, D3] = computeExpDividedDifferences(a, b, tau)
%
%   For f(x) = exp(x*tau) and every pair of nodes (a, b), returns the
%   divided differences
%       D1 = f[a, b]       = int_0^tau exp(a*(tau-u)) exp(b*u) du
%       D2 = f[a, b, 0]    = int_0^tau g(a, tau-u) exp(b*u) du
%       D3 = f[a, b, 0, 0] = int_0^tau g(a, tau-u) g(b, u) du
%   with g(lam, t) = (exp(lam*t) - 1)/lam. They are the pairwise time
%   integrals of two modes of dy/dt = A*y + c used by the adjoint and
%   tangent solves (computeAdjointGradient, computeTangentSequence), and
%   they are symmetric in a and b.
%
% INPUTS
%   a, b   arrays of nodes (eigenvalues, 1/s, possibly complex) of
%          compatible sizes (e.g. a column and a row vector: all pairs)
%   tau    scalar time (s), tau >= 0
%
% OUTPUTS
%   D1, D2, D3  arrays of the broadcast size of a and b (units s, s^2, s^3)
%
% NOTES
%   - Accurate for equal or nearly equal nodes and for nodes at or near
%     zero, where the textbook quotients cancel:
%       D1 = tau*exp(m*tau)*phi1(d*tau), m the node with the larger real
%          part, d = (other node) - m (so no overflow), phi1(z) =
%          expm1(z)/z, by its Taylor series for |z| < 1e-2;
%       D2, D3 (max(|a|,|b|)*tau >= 1): recurrences on the node p of
%          larger modulus, D2 = (D1 - g(q))/p, D3 = (D2 - h(q))/p, with q
%          the other node and g, h from computeModalIntegrals (no
%          cancellation, since |p|*tau >= 1);
%       D2, D3 (both nodes small): the series
%          f[x_0..x_m] = sum_k tau^(m+k) h_k(x)/(m+k)!, where h_k(a, b) =
%          sum_i a^i b^(k-i) (complete homogeneous polynomial; zero nodes
%          add nothing), 24 terms (truncation < 1e-24).
%   - Checked against the matrix exponential of the bidiagonal (Opitz)
%     matrix in tests/unit/AdjointTest.
%
% SEE ALSO
%   computeModalIntegrals, computeAdjointGradient, computeTangentSequence
%
% AUTHORS
%   Rishav Mallick (EOS), 2026

sz = size(a + b);
a = a + zeros(sz); b = b + zeros(sz);
if tau == 0
    D1 = zeros(sz); D2 = zeros(sz); D3 = zeros(sz);
    return
end

% ---- D1 = f[a, b], base on the node with the larger real part ----
swap = real(b) > real(a);
m = a; m(swap) = b(swap);
o = b; o(swap) = a(swap);
z = (o - m) * tau;
phi1 = ones(sz) + 0*z;
small = abs(z) < 1e-2;
zs = z(small);
ps = zeros(size(zs));
fk = 1;                                 % (k+1)!
for k = 0:7
    ps = ps + zs.^k / fk;
    fk = fk * (k + 2);
end
phi1(small) = ps;
phi1(~small) = expm1(z(~small)) ./ z(~small);
D1 = tau * exp(m * tau) .* phi1;
if nargout < 2
    return
end

% ---- D2, D3: recurrences for large nodes, series for small ones ----
big = abs(a) >= abs(b);
p = b; p(big) = a(big);                 % node of larger modulus
q = a; q(big) = b(big);
rec = abs(p) * tau >= 1;

D2 = zeros(sz) + 0*D1;
D3 = D2;
[~, gq, hq] = computeModalIntegrals(q(rec), tau);
D2(rec) = (D1(rec) - gq) ./ p(rec);
D3(rec) = (D2(rec) - hq) ./ p(rec);

% series: f[a,b,0] = sum_k tau^(k+2) h_k/(k+2)!, f[a,b,0,0] = ... /(k+3)!
as = a(~rec) * tau; bs = b(~rec) * tau;
hk = ones(size(as));                    % h_0
bk = ones(size(bs));                    % b^k
s2 = zeros(size(as)); s3 = s2;
f2 = 2; f3 = 6;                         % (k+2)!, (k+3)!
for k = 0:23
    if k > 0
        bk = bk .* bs;
        hk = as .* hk + bk;
        f2 = f2 * (k + 2); f3 = f3 * (k + 3);
    end
    s2 = s2 + hk / f2;
    s3 = s3 + hk / f3;
end
D2(~rec) = tau^2 * s2;
D3(~rec) = tau^3 * s3;
end
