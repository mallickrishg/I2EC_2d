function [V,e22dot,e23dot,parts] = extractRatesFromLinearState(sol,locked,N,M)
% EXTRACTRATESFROMLINEARSTATE  Split linear-solver states into parts.
%
%   [V, e22dot, e23dot] = extractRatesFromLinearState(sol, locked, N, M)
%   [V, e22dot, e23dot, parts] = extractRatesFromLinearState(sol, sys)
%
%   Splits states of the semi-analytical solver, ordered as
%   [creeping-patch slip rates; e22 rates; e23 rates], into fault slip
%   rate (zero on locked patches) and the two deviatoric strain-rate
%   components. It works the same way on time-integrated states (slip and
%   strain), e.g. the second output of solveAnalyticIvp.
%
%   With the system struct of assembleLinearSystem instead of (locked, N,
%   M), sol may hold Kelvin (transient) rates as well: V, e22dot and
%   e23dot are then the total slip and strain rates (Maxwell plus Kelvin),
%   and parts holds each element's share.
%
% INPUTS
%   sol     [n x Nt] states, one time per column, n = Nc + 2M (or the
%           size of sys.A)
%   locked  [N x 1] logical, true for locked fault patches
%   N       number of fault patches (rcv.N)
%   M       number of shear-zone cells (shz.N)
%   sys     struct from assembleLinearSystem (replaces locked, N, M)
%
% OUTPUTS
%   V       [N x Nt] fault slip rate (m/s); zero on locked patches
%   e22dot  [M x Nt] deviatoric strain rate e22 (1/s), total
%   e23dot  [M x Nt] shear strain rate e23 (1/s), total
%   parts   (with sys) struct: .e22M, .e23M [M x Nt] Maxwell element;
%           .e22K, .e23K [M x p x Nt] Kelvin elements (0 where a cell has
%           none); .vM [N x Nt] fault dashpot and .vK [N x pf x Nt] fault
%           Kelvin elements (0 where a patch has none, and on locked
%           patches)
%
% EXAMPLE
%   [V, e22dot, e23dot] = ...
%       extractRatesFromLinearState(sol, locked, rcv.N, shz.N);
%
% SEE ALSO
%   solveAnalyticIvp, assembleStressKernel, assembleLinearSystem
%
% AUTHORS
%   Rishav Mallick, 2025

if nargin == 2
    sys = locked;
    [V, e22dot, e23dot] = extractRatesFromLinearState(sys.P*sol, sys.locked, sys.N, sys.M);
    if nargout > 3
        Nt = size(sol, 2); M = sys.M; Nc = sys.Nc;
        parts.e22M = sol(Nc+1:Nc+M, :);
        parts.e23M = sol(Nc+M+1:Nc+2*M, :);
        parts.e22K = zeros(M, sys.p, Nt);
        parts.e23K = zeros(M, sys.p, Nt);
        for k = 1:sys.p
            kv = sys.kelvin(k);
            parts.e22K(kv.cells, k, :) = reshape(sol(kv.rows22, :), [], 1, Nt);
            parts.e23K(kv.cells, k, :) = reshape(sol(kv.rows23, :), [], 1, Nt);
        end
        parts.vM = zeros(sys.N, Nt);
        parts.vM(~sys.locked, :) = sol(1:Nc, :);
        pf = 0; if isfield(sys, 'pf'), pf = sys.pf; end
        parts.vK = zeros(sys.N, pf, Nt);
        for k = 1:pf
            fk = sys.faultKelvin(k);
            parts.vK(fk.patches, k, :) = reshape(sol(fk.rows, :), [], 1, Nt);
        end
    end
    return
end

Nt = length(sol(1,:));
V = zeros(N,Nt);
V(~locked,:) = sol(1:length(find(~locked)),:);
e22dot = sol(length(find(~locked))+1:length(find(~locked))+M,:);
e23dot = sol(length(find(~locked))+M+1:end,:);

end
