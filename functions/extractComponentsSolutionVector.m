function [V,e22dot,e23dot] = extractComponentsSolutionVector(sol,locked,N,M)
% EXTRACTCOMPONENTSSOLUTIONVECTOR  Extract velocity and strain-rate components from solution vector
%
%   [V, e22dot, e23dot] = EXTRACTCOMPONENTSSOLUTIONVECTOR(sol, locked, N, M)
%
%   Given a solution matrix `sol` from a PDE/ODE solver, this function
%   separates the solution into velocity and strain-rate components.
%
%   Inputs:
%     sol    - (n x Nt) solution matrix where columns correspond to time steps
%     locked - (N x 1) logical vector indicating which velocity degrees of freedom are constrained
%     N      - total number of velocity DOFs
%     M      - number of elements for strain-rate components
%
%   Outputs:
%     V       - (N x Nt) velocity matrix; locked DOFs remain zero
%     e22dot  - (M x Nt) first strain-rate component
%     e23dot  - (M x Nt) second strain-rate component
%
%   Notes:
%     - Assumes `sol` is ordered as [unlocked velocities; e22dot; e23dot] at each time step.
%     - The function automatically maps the unlocked DOFs back into the full velocity vector.
%
%   Example:
%       [V, e22dot, e23dot] = extractComponentsSolutionVector(sol, locked, N, M);
%
%   See also logical, zeros

Nt = length(sol(1,:));
V = zeros(N,Nt);
V(~locked,:) = sol(1:length(find(~locked)),:);
e22dot = sol(length(find(~locked))+1:length(find(~locked))+M,:);
e23dot = sol(length(find(~locked))+M+1:end,:);

end