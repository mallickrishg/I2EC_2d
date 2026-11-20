function cyc = newtonianCycle(rcv, shz, evl, Vpl, yr)
% NEWTONIANCYCLE  Rheology and periodic solution of the Newtonian example.
%
%   cyc = newtonianCycle(rcv, shz, evl, Vpl, yr)
%
%   Sets the rheology of examples 2-4 on rcv and shz and solves the
%   periodic earthquake cycle with the semi-analytical solver:
%     - mantle: Newtonian (n = 1), uniform viscosity 1e19 Pa s
%     - fault: linear viscous, tau = eta' * v, with eta' = 5e17 Pa s
%       divided by the fault width (rcv.Asigma, rcv.faultLaw = 'linear',
%       so the numerical solver integrates the same fault law)
%     - locked above 30 km depth; one earthquake every 500 years at t = 0
%       with the slip deficit of a cycle (500 yr x Vpl) on the locked zone
%     - long-term loading: plate rate on the fault, corner-flow strain
%       rates in the mantle (getStrainratesLongterm)
%
% INPUTS
%   rcv, shz  geometry.receiver (megathrust) and geometry.shearZoneReceiver
%             (mantle); handle objects, modified in place
%   evl       stress kernels (computeAllStressKernelsBem, default options)
%   Vpl       plate rate (m/s)
%   yr        seconds per year
%
% OUTPUTS
%   cyc  struct with fields
%        .Trecur   cycle length (s)
%        .locked   [N x 1] logical, locked patches
%        .slip     [N x 1] coseismic slip (m)
%        .A, .c    operator and forcing of dy/dt = A*y + c, with
%                  y = [v on creeping patches; e22dot; e23dot]
%        .b        coseismic jump in y
%        .y0       periodic state just after the earthquake
%        .eta      viscosity vector of y (MPa s/m for the fault, MPa s)
%        .solve    @(tvec) struct of V [N x nt] (0 on locked patches),
%                  e22dot, e23dot [M x nt], slip, e22, e23 (integrated
%                  from t = 0, relative to t = 0) at times tvec after the
%                  earthquake
%        .stress_change  the earthquake as a stress_change struct for
%                  viscoeq_sequence_simulation / runImposedEarthquakeCycles
%                  (Nevents = 1, Timing = 0)
%
% SEE ALSO
%   exampleModel, computePeriodicInitialConditionBalanced,
%   analyticSolveIVP, assembleStressKernel

cyc.Trecur = 500*yr;

% fault: linear viscous, eta' = 5e17 Pa s / fault width (MPa s/m)
rcv.faultLaw = 'linear';
rcv.Asigma = 1e-6*(5e17/sum(rcv.W)) .* ones(rcv.N, 1);
cyc.locked = abs(rcv.xc(:,2)) < 30e3;
rcv.pinnedPosition = cyc.locked;
rcv.Vpl = Vpl .* ones(rcv.N, 1);

% mantle: Newtonian, 1e19 Pa s everywhere
shz.n = ones(shz.N, 1);
shz.alpha = 1/(1e19*1e-6) .* ones(shz.N, 1);
[e22_lt, e23_lt] = getStrainratesLongterm(shz, rcv.dip(1)*pi/180, [0, 20e3], [-140e3, 35e3], 1);
shz.e22pl = e22_lt .* Vpl;
shz.e23pl = -e23_lt .* Vpl;

% earthquake: the slip deficit of one cycle on the locked zone
locked = cyc.locked;
cyc.slip = cyc.Trecur*Vpl .* locked;
cyc.stress_change = struct('Nevents', 1, 'Timing', 0, ...
    'dtau', evl.KK*cyc.slip, ...
    'dsigma22', evl.KL(:,:,1)*cyc.slip, 'dsigma23', evl.KL(:,:,2)*cyc.slip);

% linear system dy/dt = A*y + c, y = [v(creeping); e22dot; e23dot]
cyc.eta = [rcv.Asigma(~locked); 1./shz.alpha; 1./shz.alpha];
cyc.A = assembleStressKernel(evl, locked) ./ cyc.eta;
vlock = Vpl .* locked;
c_locking = -[evl.KK(~locked,:)*vlock; evl.KL(:,:,1)*vlock; evl.KL(:,:,2)*vlock] ./ cyc.eta;
cyc.c = c_locking - cyc.A*[rcv.Vpl(~locked); shz.e22pl; shz.e23pl];
sc = cyc.stress_change;
cyc.b = [sc.dtau(~locked); sc.dsigma22; sc.dsigma23] ./ cyc.eta;

% periodic cycle: y(T) + b = y(0)
cyc.y0 = computePeriodicInitialConditionBalanced(cyc.A, cyc.c, cyc.b, cyc.Trecur);
cyc.solve = @(tvec) evaluate(cyc, locked, rcv.N, shz.N, tvec);
end

function out = evaluate(cyc, locked, N, M, tvec)
[Y, U] = analyticSolveIVP(cyc.A, cyc.c, cyc.y0, tvec);
[out.V, out.e22dot, out.e23dot] = extractComponentsSolutionVector(Y, locked, N, M);
[out.slip, out.e22, out.e23] = extractComponentsSolutionVector(U, locked, N, M);
end
