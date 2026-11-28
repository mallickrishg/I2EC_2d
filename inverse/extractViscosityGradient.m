function gp = extractViscosityGradient(g, sys, varargin)
% EXTRACTVISCOSITYGRADIENT  Map row gradients to fault and mantle parameters.
%
%   gp = extractViscosityGradient(g, sys)
%   gp = extractViscosityGradient(g, sys, 'Parameters', 'raw')
%
%   Sums the per-row gradients of computeAdjointGradient into gradients
%   with respect to the physical parameters of each fault patch and
%   mantle cell. Maxwell: log eta_M per cell. Burgers (Kelvin elements,
%   sys.p >= 1) has three per cell and Kelvin column, which by default are
%       log eta_M,  log tauK = log(etaK/muK),  log(muK/G)
%   (relaxation time and relative strength of the Kelvin element: the
%   combinations surface data resolve), or with 'raw' log eta_M, log etaK,
%   log muK. Both are linear maps of the raw logs, so
%       dJ/dlog tauK = gEtaK,   dJ/dlog(muK/G) = gEtaK + gMuK.
%   G is any fixed modulus; it does not enter the gradient. To hold muK
%   fixed and invert eta_M and eta_K, use .etaM and .tauK and ignore
%   .muKrel (with muK fixed, a change of log tauK is the same change of
%   log etaK).
%   A Burgers fault (fault Kelvin elements, sys.pf >= 1) has the same
%   three per creeping patch and column: log eta' (.fault), and
%   log tauK = log(etaK/muK), log(muK/k0) of the fault element (k0 any
%   fixed stiffness), or with 'raw' log etaK, log muK.
%   assembleRowPerturbation is the transpose map (directions).
%
% INPUTS
%   g      struct from computeAdjointGradient (.eta, .muK) [n x r]
%   sys    struct from assembleLinearSystem
%   'Parameters'  'default' (log eta_M, log tauK, log muK/G) or 'raw'
%                 (log eta_M, log etaK, log muK)  ['default']
%
% OUTPUTS
%   gp     struct with fields ([.. x r] for r columns of g)
%          .fault   [N x r] dJ/dlog(eta') per fault patch (0 on locked)
%          .etaM    [M x r] dJ/dlog(eta_M) per cell
%          .tauK, .muKrel  [M x p x r] ('default', Kelvin only)
%          .etaK, .muK     [M x p x r] ('raw', Kelvin only)
%          .faultTauK, .faultMuKrel  [N x pf x r] ('default', fault
%                   Kelvin only)
%          .faultEtaK, .faultMuK     [N x pf x r] ('raw', fault Kelvin
%                   only)
%          Kelvin entries are 0 for cells (patches) without that element,
%          and on locked patches.
%
% SEE ALSO
%   computeAdjointGradient, assembleRowPerturbation, assembleLinearSystem
%
% AUTHORS
%   Rishav Mallick (EOS), 2026

kind = parseKind(varargin{:});
r = size(g.eta, 2);
N = sys.N; M = sys.M; Nc = sys.Nc; p = sys.p;
gp.fault = zeros(N, r);
gp.fault(~sys.locked, :) = g.eta(1:Nc, :);
gp.etaM = g.eta(Nc + (1:M), :) + g.eta(Nc + M + (1:M), :);
pf = 0; if isfield(sys, 'pf'), pf = sys.pf; end
if pf > 0
    gEtaKf = zeros(N, pf, r); gMuKf = zeros(N, pf, r);
    for k = 1:pf
        fk = sys.faultKelvin(k);
        gEtaKf(fk.patches, k, :) = reshape(g.eta(fk.rows, :), [], 1, r);
        gMuKf(fk.patches, k, :) = reshape(g.muK(fk.rows, :), [], 1, r);
    end
    if strcmp(kind, 'raw')
        gp.faultEtaK = gEtaKf; gp.faultMuK = gMuKf;
    else
        gp.faultTauK = gEtaKf; gp.faultMuKrel = gEtaKf + gMuKf;
    end
end
if p == 0
    return
end
gEtaK = zeros(M, p, r); gMuK = zeros(M, p, r);
for k = 1:p
    kv = sys.kelvin(k);
    gEtaK(kv.cells, k, :) = reshape(g.eta(kv.rows22, :) + g.eta(kv.rows23, :), [], 1, r);
    gMuK(kv.cells, k, :) = reshape(g.muK(kv.rows22, :) + g.muK(kv.rows23, :), [], 1, r);
end
if strcmp(kind, 'raw')
    gp.etaK = gEtaK; gp.muK = gMuK;
else
    gp.tauK = gEtaK; gp.muKrel = gEtaK + gMuK;
end
end

function kind = parseKind(varargin)
ip = inputParser;
ip.addParameter('Parameters', 'default', @(s) any(strcmpi(s, {'default', 'raw'})));
ip.parse(varargin{:});
kind = lower(char(ip.Results.Parameters));
end
