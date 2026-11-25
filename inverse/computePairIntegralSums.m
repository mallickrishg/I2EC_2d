function S = computePairIntegralSums(lam, tau, terms, rho)
% COMPUTEPAIRINTEGRALSUMS  Sums over time segments of pairwise modal integrals, at any time sampling.
%
%   S = computePairIntegralSums(lam, tau, terms)
%   S = computePairIntegralSums(lam, tau, terms, rho)
%
%   For eigenvalues lam and segment lengths tau(s), returns the n x n
%   matrix
%       S(i,j) = sum over terms t, segments s of
%                X_t(i,s) * Y_t(j,s) * D_{k_t}(lam_i, lam_j, tau(s))
%   with D_1, D_2, D_3 the divided differences of exp(x*tau) on the nodes
%   (lam_i, lam_j), (lam_i, lam_j, 0) and (lam_i, lam_j, 0, 0), i.e. the
%   time integrals of two modes over a segment
%   (computeExpDividedDifferences). This is the sum that
%   computeAdjointGradient forms; its cost does not depend on how many
%   distinct segment lengths there are.
%
% INPUTS
%   lam    [n x 1] eigenvalues (1/s), possibly complex, Re <= 0
%   tau    [1 x ns] segment lengths (s), >= 0
%   terms  struct array with fields
%          .k   order 1, 2 or 3
%          .X   [n x ns] coefficient of mode i in segment s
%          .Y   [n x ns] coefficient of mode j in segment s
%   rho    relative separation below which a pair of eigenvalues is
%          treated as close  [1e-4]
%
% OUTPUTS
%   S      [n x n] the sum
%
% NOTES
%   - Separable form. With F_k(x) the k-th modal integral minus its value
%     at x = 0 (shiftedModalIntegrals: expm1(x*tau), g - tau,
%     h - tau^2/2), D_k(lam_i, lam_j, tau) = (F_k(lam_j) - F_k(lam_i)) /
%     (lam_j - lam_i). The sum over segments then factors into two matrix
%     products, S = C .* (X*(Y.*F).' - (X.*F)*Y.'), C = 1./(lam_j - lam_i),
%     so all segments cost two n x n x ns products, whatever their
%     lengths. Subtracting the value at x = 0 removes the cancellation
%     between nearly equal terms for slow modes: the rounding error is
%     about eps*max(|lam_i|, |lam_j|)/|lam_i - lam_j| of the natural
%     scale tau^k.
%   - Close pairs (|lam_i - lam_j| <= rho*max(|lam_i|, |lam_j|), the
%     diagonal included) would lose digits that way; they are summed
%     directly with computeExpDividedDifferences, segment by segment,
%     which costs O(number of close pairs) per segment. With rho = 1e-4
%     the separable part keeps ~1e-12 of the natural scale.
%   - Checked against the direct per-segment sum in
%     tests/unit/AdjointTest, for spread, clustered, repeated and complex
%     spectra and log-spaced segments.
%
% SEE ALSO
%   computeExpDividedDifferences, computeAdjointGradient,
%   computeTangentSequence
%
% AUTHORS
%   Rishav Mallick (EOS), 2026

if nargin < 4 || isempty(rho)
    rho = 1e-4;
end
lam = lam(:); tau = tau(:).';
[C, I, J] = separablePairs(lam, rho);
F = cell(1, 3);
S = zeros(numel(lam)) + 0*C;
for t = 1:numel(terms)
    k = terms(t).k;
    if isempty(F{k})
        F{k} = shiftedModalIntegrals(lam, tau, k);
    end
    X = terms(t).X; Y = terms(t).Y;
    S = S + (X * (Y .* F{k}).' - (X .* F{k}) * Y.');
end
S = C .* S;

% close pairs: direct, segment by segment
if ~isempty(I)
    idx = sub2ind(size(S), I, J);
    acc = zeros(numel(I), 1) + 0*S(1);
    for s = 1:numel(tau)
        if tau(s) == 0, continue, end
        D = cell(1, 3);
        [D{1}, D{2}, D{3}] = computeExpDividedDifferences(lam(I), lam(J), tau(s));
        for t = 1:numel(terms)
            acc = acc + terms(t).X(I, s) .* terms(t).Y(J, s) .* D{terms(t).k};
        end
    end
    S(idx) = S(idx) + acc;
end
end
