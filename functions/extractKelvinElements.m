function [etaK, muK, p] = extractKelvinElements(shz)
% EXTRACTKELVINELEMENTS  Validated Kelvin elements (transient creep) of a mantle mesh.
%
%   [etaK, muK, p] = extractKelvinElements(shz)
%
%   Returns the Kelvin elements that each mantle cell has in series with its
%   steady-state (Maxwell, power-law) dashpot. Each element obeys
%       sigma = muK * epsK + etaK * depsK/dt
%   in both deviatoric components, and its strain rate adds to that of the
%   Maxwell element. One element per cell is a Burgers body; several form a
%   Prony series, e.g. for Andrade creep (computeAndradePronySeries).
%
% INPUTS
%   shz     geometry.shearZoneReceiver object (or struct) with .N and,
%           optionally, .etaK and .muK, [M x p] or [1 x p] (expanded to all
%           cells). Missing or empty properties (objects saved before they
%           existed) mean no Kelvin elements.
%
% OUTPUTS
%   etaK    [M x p] Kelvin viscosities (MPa s); Inf where a cell has no
%           element in that column
%   muK     [M x p] Kelvin moduli (MPa); set to 0 where etaK is Inf
%   p       number of columns (0: Maxwell everywhere)
%
% NOTES
%   - Errors ('extractKelvinElements:invalid') if the sizes differ, if a
%     viscosity is not positive, or if a modulus is not positive and finite
%     where the viscosity is finite (muK = 0 would be a second dashpot in
%     series: a neutral mode with no relaxed state).
%   - The viscosity convention is that of shz.alpha (edot = alpha*sigma,
%     alpha = 1/eta for n = 1), so a Kelvin element and the Maxwell dashpot
%     with the same eta have the same strain rate under the same stress.
%
% SEE ALSO
%   assembleLinearSystem, odeViscoelastic, computeAndradePronySeries,
%   geometry.shearZoneReceiver

M = shz.N;
etaK = readField(shz, 'etaK');
muK = readField(shz, 'muK');
if isempty(etaK) && isempty(muK)
    etaK = zeros(M, 0); muK = zeros(M, 0); p = 0;
    return
end
if size(etaK, 1) == 1, etaK = repmat(etaK, M, 1); end
if size(muK, 1) == 1,  muK = repmat(muK, M, 1);   end
if ~isequal(size(etaK), size(muK)) || size(etaK, 1) ~= M
    error('extractKelvinElements:invalid', ...
        'shz.etaK and shz.muK must both be [M x p] (or [1 x p]) with M = %d.', M);
end
active = isfinite(etaK);
if any(etaK(:) <= 0 | isnan(etaK(:)))
    error('extractKelvinElements:invalid', 'shz.etaK must be positive (Inf: no element).');
end
if any(~(muK(active) > 0 & isfinite(muK(active))))
    error('extractKelvinElements:invalid', ...
        'shz.muK must be positive and finite wherever shz.etaK is finite.');
end
muK(~active) = 0;
p = size(etaK, 2);
end

function v = readField(shz, name)
v = [];
if (isobject(shz) && isprop(shz, name)) || (isstruct(shz) && isfield(shz, name))
    v = shz.(name);
end
end
