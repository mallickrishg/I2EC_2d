function [dTheta, dPhi] = assembleRowPerturbation(sys, dp, varargin)
% ASSEMBLEROWPERTURBATION  Map fault and mantle parameter directions to state rows.
%
%   [dTheta, dPhi] = assembleRowPerturbation(sys, dp)
%   [dTheta, dPhi] = assembleRowPerturbation(sys, dp, 'Parameters', 'raw')
%
%   Transpose of extractViscosityGradient: turns perturbations of the
%   physical parameters (log eta' per fault patch, log eta_M per cell,
%   and for Kelvin elements log tauK and log(muK/G), or log etaK and
%   log muK with 'raw') into the per-row perturbations dTheta of
%   log(sys.eta) and dPhi of log(muK) taken by computeTangentSequence.
%   For any g: sum of gp.(f).*dp.(f) over the fields equals
%   g.eta'*dTheta + g.muK'*dPhi.
%
%   To hold a parameter fixed, leave its field out (or zero): it is then
%   not perturbed. With the default parameters, leaving out .muKrel holds
%   muK fixed while .tauK still moves etaK (etaK = tauK*muK), which is how
%   to invert eta_M and eta_K with muK known; inverting raw etaK and muK
%   instead trades them off almost perfectly (example 6).
%
% INPUTS
%   sys    struct from assembleLinearSystem
%   dp     struct with any of the fields (missing fields are zero; r
%          directions in the last dimension):
%          .fault [N x r], .etaM [M x r],
%          'default': .tauK, .muKrel [M x p x r]
%          'raw':     .etaK, .muK    [M x p x r]
%   'Parameters'  'default' or 'raw'  ['default']
%
% OUTPUTS
%   dTheta  [n x r] perturbation of log(eta) per state row
%   dPhi    [n x r] perturbation of log(muK) per state row (Kelvin rows)
%
% SEE ALSO
%   extractViscosityGradient, computeTangentSequence
%
% AUTHORS
%   Rishav Mallick (EOS), 2026

ip = inputParser;
ip.addParameter('Parameters', 'default', @(s) any(strcmpi(s, {'default', 'raw'})));
ip.parse(varargin{:});
raw = strcmpi(ip.Results.Parameters, 'raw');

N = sys.N; M = sys.M; Nc = sys.Nc; p = sys.p;
n = numel(sys.eta);
f = fieldnames(dp);
r = 1;
for i = 1:numel(f)
    v = dp.(f{i});
    if ismember(f{i}, {'fault', 'etaM'}), r = max(r, size(v, 2)); else, r = max(r, size(v, 3)); end
end
get2 = @(name, rows) fieldOr(dp, name, zeros(rows, r));
get3 = @(name) fieldOr(dp, name, zeros(M, p, r));

dTheta = zeros(n, r); dPhi = zeros(n, r);
df = get2('fault', N);
dTheta(1:Nc, :) = df(~sys.locked, :);
dm = get2('etaM', M);
dTheta(Nc + (1:M), :) = dm;
dTheta(Nc + M + (1:M), :) = dm;
if p == 0
    return
end
if raw
    dEtaK = get3('etaK'); dMuK = get3('muK');
else
    dTau = get3('tauK'); dRel = get3('muKrel');
    dEtaK = dTau + dRel; dMuK = dRel;     % etaK = tauK*muK, muK = G*muKrel
end
for k = 1:p
    kv = sys.kelvin(k);
    de = reshape(dEtaK(kv.cells, k, :), [], r);
    dmu = reshape(dMuK(kv.cells, k, :), [], r);
    dTheta(kv.rows22, :) = de; dTheta(kv.rows23, :) = de;
    dPhi(kv.rows22, :) = dmu; dPhi(kv.rows23, :) = dmu;
end
end

function v = fieldOr(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end
