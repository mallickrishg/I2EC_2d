function sys = assembleLinearSystem(evl, rcv, shz, locked, stress_change)
% ASSEMBLELINEARSYSTEM  Linear cycle ODE dy/dt = A*y + c for Maxwell, Burgers or Andrade mantles.
%
%   sys = assembleLinearSystem(evl, rcv, shz, locked)
%   sys = assembleLinearSystem(evl, rcv, shz, locked, stress_change)
%
%   Builds the operator A, the forcing c and the coseismic jump b of the
%   semi-analytical solver for a linear viscous fault and a linear (n = 1)
%   mantle whose cells may have Kelvin elements in series with the Maxwell
%   dashpot (shz.etaK, shz.muK; see extractKelvinElements). A, c and b go
%   unchanged into computePeriodicInitialCondition and solveAnalyticIvp.
%
%   The state is made of rates only:
%       y = [v on creeping patches; edotM22; edotM23;
%            edotK22, edotK23 of Kelvin column 1 on the cells that have it;
%            ... of column 2; ...]
%   Each entry r maps to one entry map(r) of the Maxwell state
%   x = [v; e22dot; e23dot] (n0 = Nc + 2M entries), and the kernels act on
%   the total rates x = P*y (Maxwell plus Kelvin):
%       dsigma/dt = K*(P*y - y_pl) + sigma_lock,
%       K = assembleStressKernel(evl, locked),
%   and each rate follows the stress of its element,
%       dy_r/dt = (dsigma/dt)_map(r) / eta_r - (muK/etaK)_r * y_r,
%   so that
%       A = K(map,map) ./ eta - diag(decay),
%       c = sigma_lock(map) ./ eta - (K(map,:) ./ eta)*y_pl,
%       b = dsigma_coseismic(map) ./ eta.
%   With no Kelvin elements map = 1:n0 and these are the expressions of the
%   Maxwell solver (examples/newtonianCycle) operation for operation.
%
% INPUTS
%   evl            stress kernels (.KK, .KL, .LK, .LL) from
%                  computeAllStressKernelsBem
%   rcv            geometry.receiver: .N, .Vpl (m/s), .Asigma (fault
%                  viscosity per unit length, MPa s/m)
%   shz            geometry.shearZoneReceiver: .N, .alpha (1/(MPa s)),
%                  .n (must be 1), .e22pl, .e23pl (1/s), and optionally
%                  .etaK, .muK
%   locked         [N x 1] logical, true for locked patches (v = 0; they
%                  load the rest at the plate rate)
%   stress_change  optional struct of events (.dtau [N x K], .dsigma22,
%                  .dsigma23 [M x K]), as for integrateImposedCycles
%
% OUTPUTS
%   sys   struct with fields
%         .A       [n x n] operator (1/s)
%         .c       [n x 1] forcing
%         .b       [n x K] jump of y at each event (empty without
%                  stress_change)
%         .eta     [n x 1] viscosity of the element of each row (MPa s/m
%                  on the fault, MPa s in the mantle)
%         .decay   [n x 1] muK/etaK on Kelvin rows, 0 elsewhere (1/s)
%         .map     [n x 1] entry of x = [v; e22dot; e23dot] of each row
%         .P       [n0 x n] sparse sum: total rates x = P*y
%         .kelvin  [1 x p] struct array, per Kelvin column: .cells (cell
%                  indices that have it), .rows22, .rows23 (their rows in y)
%         .locked, .N, .M, .Nc, .n0, .p   sizes and the locked mask
%
% NOTES
%   - Kelvin strain: epsK = (sigma - etaK*edotK)/muK with sigma = edotM/alpha.
%     It is continuous through an event, so over a periodic cycle the
%     Kelvin rates integrate to zero and the cycle-mean slip and strain
%     rates are those of the Maxwell model.
%   - Springs and dashpots form a passive network: the spectrum stays real
%     and negative when the Maxwell operator's is. Cells with etaK = Inf
%     carry no Kelvin state, so etaK -> Inf reduces to Maxwell.
%   - Size n = Nc + 2M + 2*(number of finite etaK); eig is O(n^3).
%   - Use extractRatesFromLinearState(Y, sys) to split states into slip
%     and (total, Maxwell and Kelvin) strain rates.
%
% EXAMPLE
%   shz.etaK = 1e17*1e-6*ones(shz.N,1); shz.muK = 3e4*ones(shz.N,1);  % Burgers
%   sys = assembleLinearSystem(evl, rcv, shz, locked, stress_change);
%   y0 = computePeriodicInitialCondition(sys.A, sys.c, sys.b, Trecur);
%   [V, e22dot, e23dot] = extractRatesFromLinearState( ...
%       solveAnalyticIvp(sys.A, sys.c, y0, tvec), sys);
%
% SEE ALSO
%   assembleStressKernel, extractKelvinElements, computeAndradePronySeries,
%   computePeriodicInitialCondition, solveAnalyticIvp,
%   extractRatesFromLinearState

if ~all(shz.n(:) == 1)
    error('assembleLinearSystem:nonlinear', ...
        'The semi-analytical solver needs a linear mantle (shz.n = 1 everywhere).');
end
locked = logical(locked(:));
N = rcv.N; M = shz.N; Nc = sum(~locked); n0 = Nc + 2*M;
[etaK, muK, p] = extractKelvinElements(shz);

% ---- rows: Maxwell state, then the Kelvin columns ----
map = (1:n0)';
eta = [rcv.Asigma(~locked); 1./shz.alpha; 1./shz.alpha];
decay = zeros(n0, 1);
kelvin = struct('cells', cell(1, p), 'rows22', [], 'rows23', []);
for k = 1:p
    cells = find(isfinite(etaK(:,k)));
    nk = numel(cells);
    kelvin(k).cells = cells;
    kelvin(k).rows22 = numel(map) + (1:nk)';
    kelvin(k).rows23 = numel(map) + nk + (1:nk)';
    map = [map; Nc + cells; Nc + M + cells]; %#ok<AGROW>
    eta = [eta; etaK(cells,k); etaK(cells,k)]; %#ok<AGROW>
    decay = [decay; repmat(muK(cells,k)./etaK(cells,k), 2, 1)]; %#ok<AGROW>
end
n = numel(map);

% ---- operator, forcing and jump ----
K = assembleStressKernel(evl, locked);
Kr = K(map,:) ./ eta;                       % stressing rate per row viscosity
sys.A = Kr(:,map);
if p > 0
    sys.A = sys.A - diag(decay);
end
vlock = rcv.Vpl .* locked;
slock = -[evl.KK(~locked,:)*vlock; evl.KL(:,:,1)*vlock; evl.KL(:,:,2)*vlock];
ypl = [rcv.Vpl(~locked); shz.e22pl; shz.e23pl];
sys.c = slock(map) ./ eta - Kr*ypl;
sys.b = [];
if nargin > 4 && ~isempty(stress_change)
    sc = stress_change;
    dsig = [sc.dtau(~locked,:); sc.dsigma22; sc.dsigma23];
    sys.b = dsig(map,:) ./ eta;
end

sys.eta = eta;
sys.decay = decay;
sys.map = map;
sys.P = sparse(map, (1:n)', 1, n0, n);
sys.kelvin = kelvin;
sys.locked = locked;
sys.N = N; sys.M = M; sys.Nc = Nc; sys.n0 = n0; sys.p = p;
end
