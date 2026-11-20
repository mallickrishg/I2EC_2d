function [V,e22dot,e23dot] = extractComponentsSolutionVector(sol,locked,N,M)
% EXTRACTCOMPONENTSSOLUTIONVECTOR  Split linear-solver states into parts.
%
%   [V, e22dot, e23dot] = extractComponentsSolutionVector(sol, locked, N, M)
%
%   Splits states of the semi-analytical solver, ordered as
%   [creeping-patch slip rates; e22 rates; e23 rates], into fault slip
%   rate (zero on locked patches) and the two deviatoric strain-rate
%   components. It works the same way on time-integrated states (slip and
%   strain), e.g. the second output of analyticSolveIVP.
%
% INPUTS
%   sol     [n x Nt] states, one time per column, n = Nc + 2M
%   locked  [N x 1] logical, true for locked fault patches
%   N       number of fault patches (rcv.N)
%   M       number of shear-zone cells (shz.N)
%
% OUTPUTS
%   V       [N x Nt] fault slip rate (m/s); zero on locked patches
%   e22dot  [M x Nt] deviatoric strain rate e22 (1/s)
%   e23dot  [M x Nt] shear strain rate e23 (1/s)
%
% EXAMPLE
%   [V, e22dot, e23dot] = ...
%       extractComponentsSolutionVector(sol, locked, rcv.N, shz.N);
%
% SEE ALSO
%   analyticSolveIVP, assembleStressKernel
%
% AUTHORS
%   Rishav Mallick, 2025

Nt = length(sol(1,:));
V = zeros(N,Nt);
V(~locked,:) = sol(1:length(find(~locked)),:);
e22dot = sol(length(find(~locked))+1:length(find(~locked))+M,:);
e23dot = sol(length(find(~locked))+M+1:end,:);

end