function KLmodes = computeKLModes(shz)
% COMPUTEKLMODES Compute Karhunen-Loève modes from spatial correlation matrix
%
%   KLmodes = COMPUTEKLMODES(shz) computes the Karhunen-Loève (KL) modes 
%   by constructing an exponential correlation matrix based on normalized 
%   Euclidean distances between spatial points and performing singular 
%   value decomposition.
%
%   Input:
%       shz     - Structure containing spatial discretization information
%                 Required fields:
%                 .N  : Number of spatial points (scalar)
%                 .xc : Coordinates of points (N x 2 matrix)
%                       Column 1: x-coordinates
%                       Column 2: y-coordinates
%
%   Output:
%       KLmodes - Karhunen-Loève modes (N x N matrix)
%                 Each column is a KL mode (spatial basis function)
%                 Modes are ordered by singular values (descending)
%                 Forms an orthonormal basis
%
%   Algorithm:
%       1. Computes pairwise Euclidean distances between all points
%       2. Normalizes distances to [0,1] range
%       3. Constructs exponential correlation matrix: C(i,j) = exp(-r_norm(i,j))
%       4. Performs SVD to extract spatial modes (left singular vectors)
%
%   Example:
%       shz.N = 100;
%       shz.xc = rand(100, 2);
%       KLmodes = computeKLModes(shz);
%       firstMode = KLmodes(:, 1);  % Most energetic mode
%
%   See also: SVD, EIG

r = zeros(shz.N,shz.N);
for i = 1:shz.N
    % compute distance
    r(:,i) = sqrt((shz.xc(:,1)-shz.xc(i,1)).^2 + (shz.xc(:,2)-shz.xc(i,2)).^2);
end
% normalize r so that r = 0 -> 1
rnorm = r./max(r(:));
C_r = exp(-rnorm);

if false
    % compute eigen values
    [eigenvectors,lambda] = eig(C_r);
    % sort eigenvalues
    [~, idx] = sort(diag(lambda), 'descend');
    KLmodes = eigenvectors(:, idx);
else
    % use SVD
    [U,~,~] = svd(C_r);
    KLmodes = U;
end

end