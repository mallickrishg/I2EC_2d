function r = drchrnd(a,n)
% DRCHRND  Draw random samples from a Dirichlet distribution.
%
%   r = drchrnd(a, n)
%
%   Draws n samples from the Dirichlet distribution with concentration
%   parameters a by normalising independent gamma variates.
%
% INPUTS
%   a       [1 x p] concentration parameters (positive)
%   n       number of samples
%
% OUTPUTS
%   r       [n x p] samples; each row is non-negative and sums to 1
%
% NOTES
%   - Needs the Statistics and Machine Learning Toolbox (gamrnd).
%
% SEE ALSO
%   gamrnd
%
% AUTHORS
%   Rishav Mallick, 2023

p = length(a);
r = gamrnd(repmat(a,n,1),1,n,p);
r = r ./ repmat(sum(r,2),1,p);
end