function fwd = solveLinearSequence(sys, seq, varargin)
% SOLVELINEARSEQUENCE  Periodic spin-up then an observed event sequence, in modal form.
%
%   fwd = solveLinearSequence(sys, seq)
%   fwd = solveLinearSequence(sys, seq, 'Balance', tf)
%
%   Forward model of the adjoint (inverse) path. The linear system
%   dy/dt = A*y + c of assembleLinearSystem is first in its periodic
%   cycle: one characteristic event (state jump b0) every T, with
%   y(0+) = y(T-) + b0. The observed window starts at age seq.age of that
%   cycle (observed clock t = 0), and then has its own events, each a
%   state jump, at times seq.tEvents. Returned are the states and the
%   integrated states (slip and strain) relative to a reference time,
%   at the observation times, plus everything the adjoint
%   (computeAdjointGradient) and tangent (computeTangentSequence) reuse:
%   one eigendecomposition of A, its inverse, and the modal states at
%   every breakpoint.
%
% INPUTS
%   sys    struct from assembleLinearSystem (.A, .c, .eta, .decay, ...)
%   seq    struct describing the time line (times in s):
%          .T         recurrence time of the spin-up cycle (> 0)
%          .b0        [n x 1] state jump of the characteristic event
%          .age       age of the spin-up cycle at observed t = 0, in
%                     [0, T]  [T: just before the next characteristic
%                     event]
%          .tEvents   [1 x K] observed event times (>= 0)  [none]
%          .bEvents   [n x K] their state jumps (e.g. columns of sys.b)
%          .tObs      [1 x nObs] observation times (>= 0)
%          .tRef      reference time of the integrated states (>= 0)
%                     [min(tObs)]
%   'Balance'  true or false: eigendecompose with or without balancing;
%              'auto' (default) picks the variant whose modal
%              propagation over one cycle is closer to the matrix
%              exponential, as solveAnalyticIvp does. Pass fwd.balance of
%              a first solve to skip the check in later solves.
%
% OUTPUTS
%   fwd    struct with fields
%          .U        [n x nObs] integral of y from tRef to tObs(j)
%          .Y        [n x nObs] y at tObs(j) (after any event at that time)
%          .y0       [n x 1] periodic state just after a characteristic
%                    event
%          .ystart   [n x 1] state at observed t = 0, before any event
%          .V, .Vi, .lam   modes: A = V*diag(lam)*Vi
%          .balance  the eigendecomposition variant used
%          .bp, .zs, .bhat, .ch, .z0, .za, .b0hat, .iObs, .iRef
%                    breakpoints and modal states (for the adjoint and
%                    tangent)
%          .seq, .sys   the inputs
%
% NOTES
%   - Breakpoints are 0, the event times, the observation times and tRef;
%     between them every mode evolves exactly (computeModalIntegrals).
%   - The periodic state is solved in modal form,
%     z0 = (g(T)*V\c + V\b0) ./ (-expm1(lam*T)), with expm1, so slow
%     modes keep their digits. computePeriodicInitialCondition solves the
%     same condition directly; the two agree to the eigendecomposition's
%     accuracy (tests/unit/AdjointTest).
%   - Growing and neutral modes are errors, as in
%     computePeriodicInitialCondition.
%   - Events at the same time add their jumps; an event at t = 0 acts after
%     the start state.
%
% EXAMPLE
%   sys = assembleLinearSystem(evl, rcv, shz, locked, stress_change);
%   seq = struct('T', 500*yr, 'b0', sys.b(:,1), 'tEvents', 0, ...
%                'bEvents', sys.b(:,1), 'tObs', (0.1:0.5:20)*yr);
%   fwd = solveLinearSequence(sys, seq);
%   u = H*fwd.U;           % H from assembleStationKernels
%
% SEE ALSO
%   assembleLinearSystem, computeAdjointGradient, computeTangentSequence,
%   computeDisplacementMisfit, computePeriodicInitialCondition,
%   solveAnalyticIvp
%
% AUTHORS
%   Rishav Mallick (EOS), 2026

ip = inputParser;
ip.addParameter('Balance', 'auto');
ip.parse(varargin{:});
balance = ip.Results.Balance;

A = sys.A; c = sys.c; n = size(A, 1);
seq = completeSequence(seq, n);
T = seq.T;

% ---- eigendecomposition (one variant) ----
if ischar(balance) || isstring(balance)
    [V, lam, nGrowing] = eigOperator(A);
    checkModes(lam, nGrowing, T);
    [Vb, lamb] = eigOperator(A, [], true);
    ytest = seq.b0;
    if ~any(ytest), ytest = c; end
    X = expm([A, c; zeros(1, n+1)] * T);
    yref = X(1:n, 1:n)*ytest + X(1:n, end);
    scale = abs(yref) + 1e-12*max(abs(yref)) + realmin;
    balance = propagationError(Vb, lamb, ytest, c, T, yref, scale) < ...
              propagationError(V, lam, ytest, c, T, yref, scale);
    if balance
        V = Vb; lam = lamb;
    end
else
    [V, lam, nGrowing] = eigOperator(A, [], logical(balance));
    checkModes(lam, nGrowing, T);
end
Vi = V \ eye(n);
ch = Vi * c;

% ---- periodic spin-up, in modal form ----
b0hat = Vi * seq.b0;
[~, gT] = computeModalIntegrals(lam, T);
z0 = (gT .* ch + b0hat) ./ (-expm1(lam * T));
[ea, ga] = computeModalIntegrals(lam, seq.age);
za = ea .* z0 + ga .* ch;

% ---- observed window: breakpoints ----
bp = unique([0, seq.tEvents(:)', seq.tObs(:)', seq.tRef]);
nbp = numel(bp);
bhat = zeros(n, nbp);
for k = 1:numel(seq.tEvents)
    i = find(bp == seq.tEvents(k), 1);
    bhat(:, i) = bhat(:, i) + Vi * seq.bEvents(:, k);
end
[~, iObs] = ismember(seq.tObs(:)', bp);
iRef = find(bp == seq.tRef, 1);

zs = zeros(n, nbp) + 0*za;              % state just after bp(i)
UZ = zs;                                % int_0^bp(i) z
z = za; uz = zeros(n, 1) + 0*za;
for i = 1:nbp
    if i > 1
        [e, g, h] = computeModalIntegrals(lam, bp(i) - bp(i-1));
        uz = uz + g .* z + h .* ch;
        z = e .* z + g .* ch;
    end
    z = z + bhat(:, i);
    zs(:, i) = z;
    UZ(:, i) = uz;
end

fwd.U = real(V * (UZ(:, iObs) - UZ(:, iRef)));
fwd.Y = real(V * zs(:, iObs));
fwd.y0 = real(V * z0);
fwd.ystart = real(V * za);
fwd.V = V; fwd.Vi = Vi; fwd.lam = lam;
fwd.balance = balance;
fwd.bp = bp; fwd.zs = zs; fwd.bhat = bhat; fwd.ch = ch;
fwd.z0 = z0; fwd.za = za; fwd.b0hat = b0hat;
fwd.iObs = iObs; fwd.iRef = iRef;
fwd.seq = seq; fwd.sys = sys;
end

function seq = completeSequence(seq, n)
if ~isfield(seq, 'age') || isempty(seq.age), seq.age = seq.T; end
if ~isfield(seq, 'tEvents') || isempty(seq.tEvents)
    seq.tEvents = zeros(1, 0); seq.bEvents = zeros(n, 0);
end
if ~isfield(seq, 'tRef') || isempty(seq.tRef), seq.tRef = min(seq.tObs); end
seq.tEvents = seq.tEvents(:)'; seq.tObs = seq.tObs(:)';
if ~(isscalar(seq.T) && seq.T > 0)
    error('solveLinearSequence:T', 'seq.T must be a positive scalar.');
end
if seq.age < 0 || seq.age > seq.T
    error('solveLinearSequence:age', 'seq.age must lie in [0, seq.T].');
end
if any([seq.tEvents, seq.tObs, seq.tRef] < 0)
    error('solveLinearSequence:times', 'Event, observation and reference times must be >= 0.');
end
if size(seq.b0, 1) ~= n || size(seq.bEvents, 1) ~= n || size(seq.bEvents, 2) ~= numel(seq.tEvents)
    error('solveLinearSequence:size', 'seq.b0 must be [n x 1] and seq.bEvents [n x numel(tEvents)].');
end
end

function checkModes(lam, nGrowing, T)
if nGrowing > 0
    error('I2EC:growingModes', ['solveLinearSequence: the operator has %d growing modes ' ...
        '(largest Re(lambda) = %.3g 1/s).'], nGrowing, max(real(lam)));
end
if any(abs(real(lam)) * T <= 1e-12)
    error('I2EC:neutralModes', ['solveLinearSequence: the operator has modes that neither ' ...
        'decay nor grow over a cycle, so the periodic state is not determined.']);
end
end

function err = propagationError(V, lambda, y0, c, t, yref, scale)
[e, g] = computeModalIntegrals(lambda, t);
y = real(V * (e .* (V \ y0) + g .* (V \ c)));
err = norm((y - yref) ./ scale) / sqrt(numel(y));
end
