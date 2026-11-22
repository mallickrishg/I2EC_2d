function fit = computeAndradePronySeries(beta, tmin, tmax, varargin)
% COMPUTEANDRADEPRONYSERIES  Kelvin elements (Prony series) that creep like Andrade.
%
%   fit = computeAndradePronySeries(beta, tmin, tmax)
%   fit = computeAndradePronySeries(beta, tmin, tmax, 'TermsPerDecade', 3)
%
%   Fits a sum of Kelvin elements to Andrade transient creep: under a
%   constant stress sigma, Andrade creep is
%       epsA(t) = beta * sigma * t^(1/3),
%   and p Kelvin elements in series creep as
%       epsP(t) = sigma * sum_k J_k * (1 - exp(-t/tau_k)),
%   J_k = 1/muK_k, tau_k = etaK_k/muK_k. The relaxation times tau_k are
%   log-spaced on [tmin, tmax] widened by 'Margin' decades at each end; the
%   compliances J_k >= 0 minimise the
%   relative error of epsP against epsA at 40 log-spaced times per decade
%   on [tmin, tmax] (lsqnonneg). The fit has a fixed shape, so it is done
%   once and scaled by beta: per-cell beta costs nothing. Assign the result
%   to the mantle mesh as
%       shz.etaK = fit.etaK; shz.muK = fit.muK;
%   and both solvers then add Andrade creep in series with the Maxwell
%   (steady-state) dashpot.
%
% INPUTS
%   beta    [M x 1] or scalar Andrade coefficient (1/(MPa s^(1/3))), in the
%           strain-rate convention of shz.alpha (edot = alpha*sigma);
%           beta = 0 gives no Kelvin elements in that cell
%   tmin    shortest time to fit (s), e.g. 1e-3 yr
%   tmax    longest time to fit (s), e.g. a few cycle lengths
%
% OPTIONS (name-value)
%   'TermsPerDecade'  relaxation times per decade  [1]
%   'Margin'          decades by which the relaxation times extend beyond
%                     [tmin, tmax]; without it the ends of the range are
%                     fitted poorly (10% error at tmin)  [0.5]
%
% OUTPUTS
%   fit   struct with fields
%         .etaK, .muK   [M x p] Kelvin viscosities (MPa s) and moduli
%                       (MPa); Inf/0 where beta = 0 or a J_k is zero
%         .tau          [1 x p] relaxation times (s)
%         .J            [1 x p] compliances per unit beta (s^(1/3)):
%                       muK = 1./(beta*J)
%         .maxRelErr    max |epsP/epsA - 1| on [tmin, tmax]
%         .rmsRelErr    rms of epsP/epsA - 1 on [tmin, tmax] (log-uniform)
%         .t, .relErr   the fit times (s) and relative errors there
%
% NOTES
%   - Outside [tmin, tmax] the series departs from Andrade by construction:
%     below tmin it creeps too slowly (t^(1/3) has an infinite initial
%     rate) and above tmax it saturates at sigma*beta*sum(J) while Andrade
%     creep grows without bound. The saturation gives a relaxed state, so
%     a periodic cycle exists; choose tmax of a few cycle lengths.
%   - Maximum relative error on [tmin, tmax] (see the table in
%     docs/dev/transient-rheology.md): about 1% with the defaults, 0.4% with
%     2 terms per decade, 1e-4 with 2 terms per decade and a margin of 1.
%     The shortest relaxation time is tmin/10^Margin: short times make the
%     numerical solver stiff (use @ode15s), and every term adds 2M states.
%   - Elements whose fitted J_k is zero are dropped (lsqnonneg is sparse).
%
% EXAMPLE
%   yr = 3.15e7;
%   fit = computeAndradePronySeries(beta*ones(shz.N,1), 1e-3*yr, 2500*yr);
%   shz.etaK = fit.etaK; shz.muK = fit.muK;
%
% SEE ALSO
%   extractKelvinElements, assembleLinearSystem, odeViscoelastic, lsqnonneg

ip = inputParser;
ip.addParameter('TermsPerDecade', 1, @(x) isscalar(x) && x > 0);
ip.addParameter('Margin', 0.5, @(x) isscalar(x) && x >= 0);
ip.parse(varargin{:});
if ~(tmin > 0 && tmax > tmin)
    error('computeAndradePronySeries:range', 'Need 0 < tmin < tmax.');
end
beta = beta(:);
if any(beta < 0 | ~isfinite(beta))
    error('computeAndradePronySeries:beta', 'beta must be finite and non-negative.');
end

decades = log10(tmax/tmin);
a = ip.Results.Margin;
p = max(2, ceil(ip.Results.TermsPerDecade*(decades + 2*a)) + 1);
tau = logspace(log10(tmin) - a, log10(tmax) + a, p);
t = logspace(log10(tmin), log10(tmax), max(50, ceil(40*decades)))';

% relative-error least squares: rows divided by t^(1/3)
G = (1 - exp(-t./tau)) ./ t.^(1/3);
J = lsqnonneg(G, ones(size(t)));
keep = J' > 0;
tau = tau(keep); J = J(keep)';

relErr = (1 - exp(-t./tau))*J' ./ t.^(1/3) - 1;
fit.tau = tau;
fit.J = J;
fit.t = t;
fit.relErr = relErr;
fit.maxRelErr = max(abs(relErr));
fit.rmsRelErr = sqrt(mean(relErr.^2));

% per cell: muK = 1/(beta*J), etaK = tau*muK
on = beta > 0;
fit.muK = zeros(numel(beta), numel(J));
fit.etaK = Inf(numel(beta), numel(J));
fit.muK(on,:) = 1 ./ (beta(on) * J);
fit.etaK(on,:) = fit.muK(on,:) .* tau;
end
