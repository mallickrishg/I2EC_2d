function d = pairDistance(p, q)
%PAIRDISTANCE  Euclidean distance between every row of P and every row of Q.
%
%   d = pairDistance(p, q)    p: [n x 2], q: [m x 2]  ->  d: [n x m]
%
% Plain MATLAB replacement for pdist2 (no Statistics Toolbox needed).

d = sqrt((p(:,1) - q(:,1)').^2 + (p(:,2) - q(:,2)').^2);

end
