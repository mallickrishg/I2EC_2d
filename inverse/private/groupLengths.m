function [ut, grp] = groupLengths(tau)
% GROUPLENGTHS  Group segment lengths that agree to 1e-12 (relative).
%
%   [ut, grp] = groupLengths(tau)
%
%   ut(grp(i)) is the representative length of tau(i) (the first of its
%   group in sorted order). Uniform observation epochs give one group
%   even when their differences differ in the last bits.
[ts, order] = sort(tau(:)');
newGroup = [true, diff(ts) > 1e-12 * ts(2:end)];
gs = cumsum(newGroup);
ut = ts(newGroup);
grp = zeros(size(ts));
grp(order) = gs;
end
