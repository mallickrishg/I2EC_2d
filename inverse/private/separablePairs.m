function [C, I, J] = separablePairs(lam, rho)
% SEPARABLEPAIRS  Cauchy matrix 1./(lam_j - lam_i) off the close pairs, and the close pairs.
%
%   [C, I, J] = separablePairs(lam, rho)
%
%   Pairs with |lam_i - lam_j| <= rho*max(|lam_i|, |lam_j|) (the diagonal
%   included) are close: C is 0 there and (I, J) lists them, to be summed
%   directly with computeExpDividedDifferences. Used by
%   computePairIntegralSums and computeTangentSequence.
lam = lam(:);
d = lam.' - lam;                                  % d(i,j) = lam_j - lam_i
close = abs(d) <= rho * max(abs(lam), abs(lam.'));
C = 1 ./ d;
C(close) = 0;
[I, J] = find(close);
end
