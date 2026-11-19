function r = reciprocityRatio(M, Mt, mask)
%RECIPROCITYRATIO  Element-wise ratio M ./ Mt over MASK (Betti check).
%
%   r = reciprocityRatio(M, Mt, mask)
%
% For a weighted kernel pair that should satisfy Betti's reciprocal
% theorem, M(mask) ./ Mt(mask) is 1.  Collocation at element centres makes
% near-field pairs inexact, so pass a mask of well-separated pairs and
% look at the median and percentiles of r.

r = M(mask) ./ Mt(mask);
r = r(isfinite(r));

end
