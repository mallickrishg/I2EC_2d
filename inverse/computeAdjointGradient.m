function g = computeAdjointGradient(fwd, S)
% COMPUTEADJOINTGRADIENT  Adjoint gradient of a misfit with respect to log viscosity per state row.
%
%   g = computeAdjointGradient(fwd, S)
%
%   Vector-Jacobian product of the forward model of solveLinearSequence.
%   For a misfit J whose derivative with respect to the integrated states
%   fwd.U is S (S(:,j) = dJ/dU(:,j)), returns dJ/dlog(eta_r) for the
%   viscosity of every state row r, and dJ/dlog(muK_r) for the Kelvin
%   spring of every Kelvin row. It knows nothing about stations: S comes
%   from computeDisplacementMisfit (S = H'*W*r) or any other misfit of U.
%   Map rows to physical parameters with extractViscosityGradient.
%
% INPUTS
%   fwd   struct from solveLinearSequence
%   S     [n x nObs] dJ/dU, one column per observation time
%
% OUTPUTS
%   g     struct with fields
%         .eta  [n x 1] dJ/dlog(sys.eta(r)), the viscosity of the
%               element of row r (fault patch, Maxwell or Kelvin dashpot)
%         .muK  [n x 1] dJ/dlog(muK) on Kelvin rows, 0 elsewhere
%
% NOTES
%   - Every row of A, c and b is divided by the viscosity of its element
%     (assembleLinearSystem), so d(A*y + c)/dlog(eta_r) = -e_r*ydot_r and
%     db/dlog(eta_r) = -e_r*b_r. With the adjoint lambda
%     (-dlambda/dt = A'*lambda + s, s(t) = sum of S(:,j) over the
%     observations whose window [tRef, tObs(j)] contains t),
%         dJ/dlog(eta_r) = -int lambda_r dy_r
%     as a Stieltjes integral: the jumps of y at events count. The spin-up
%     cycle adds the same integral with the adjoint of the periodicity
%     condition, nu = (I - expm(A*T))^(-T) * lambda_p(0), diagonal in modal
%     form. The Kelvin spring enters A only on the diagonal, so
%         dJ/dlog(muK_r) = -(muK/etaK)_r * int lambda_r y_r dt.
%   - In modal form (lambda = Vi.'*alpha, y = V*z) the time integrals are
%     pairwise divided differences (computeExpDividedDifferences). They
%     are summed into one n x n matrix Q, and dJ/dlog(eta) =
%     -rowsum((Vi.'*Q) .* V): the cost is one n x n product on top of the
%     forward's eigendecomposition, plus O(n^2) per distinct segment
%     length.
%   - Segment lengths equal to 1e-12 (relative) are grouped and share their
%     divided differences (uniform epochs: one group).
%   - Rows with a non-finite sys.eta (no viscosity) get 0.
%   - Checked against central finite differences and, through the
%     dot-product test, against computeTangentSequence
%     (tests/unit/AdjointTest).
%
% EXAMPLE
%   fwd = solveLinearSequence(sys, seq);
%   [J, res] = computeDisplacementMisfit(obs, fwd);
%   g = computeAdjointGradient(fwd, res.S);
%   gp = extractViscosityGradient(g, sys);       % per cell: gp.etaM, ...
%
% SEE ALSO
%   solveLinearSequence, computeTangentSequence, extractViscosityGradient,
%   computeDisplacementMisfit, computeExpDividedDifferences
%
% AUTHORS
%   Rishav Mallick (EOS), 2026

V = fwd.V; Vi = fwd.Vi; lam = fwd.lam; ch = fwd.ch;
n = numel(lam);
bp = fwd.bp; nbp = numel(bp); nseg = nbp - 1;
seq = fwd.seq; sys = fwd.sys;
decay = sys.decay;
doMu = any(decay ~= 0);

% ---- adjoint sources per segment ----
tObs = seq.tObs; tRef = seq.tRef;
lo = min(tObs, tRef); hi = max(tObs, tRef);
sgn = sign(tObs - tRef);
segIn = (bp(1:end-1)' >= lo) & (bp(2:end)' <= hi);   % [nseg x nObs]
shat = V.' * (S * (segIn .* sgn).');                % [n x nseg]

% ---- backward sweep: alpha at the end and start of each segment ----
aend = zeros(n, nseg) + 0*shat;
astart = aend;
alpha = zeros(n, 1) + 0*lam;
tau = diff(bp);
for i = nseg:-1:1
    aend(:, i) = alpha;
    [e, gg] = computeModalIntegrals(lam, tau(i));
    alpha = e .* alpha + gg .* shat(:, i);
    astart(:, i) = alpha;
end
mu = alpha;                                         % at observed t = 0

% ---- spin-up cycle: [0, age] to the start state, and the periodicity ----
T = seq.T; a = seq.age;
nu = exp(lam * a) .* mu ./ (-expm1(lam * T));

% ---- Q = sum of pairwise integrals; jump terms first ----
Q = astart * fwd.bhat(:, 1:nseg).' + nu * fwd.b0hat.';
Qy = zeros(n) + 0*Q;

% all segments: [observed..., spin-up to age a, whole cycle]
segTau = [tau, a, T];
segA = [aend, mu, nu];
segS = [shat, zeros(n, 2)];
segZ = [fwd.zs(:, 1:nseg), fwd.z0, fwd.z0];
[ut, grp] = groupLengths(segTau);
for u = 1:numel(ut)
    k = grp == u;
    if ut(u) == 0, continue, end
    [E1, E2, E3] = computeExpDividedDifferences(lam, lam.', ut(u));
    W = lam .* segZ(:, k) + ch;                     % zdot(t0 + s) = exp(lam*s).*w
    Q = Q + E1 .* (segA(:, k) * W.') + E2 .* (segS(:, k) * W.');
    if doMu
        Qy = Qy + E1 .* (segA(:, k) * segZ(:, k).') ...
                + E2 .* (sum(segA(:, k), 2) * ch.' + segS(:, k) * segZ(:, k).') ...
                + E3 .* (sum(segS(:, k), 2) * ch.');
    end
end

g.eta = -real(sum((Vi.' * Q) .* V, 2));
g.muK = zeros(n, 1);
if doMu
    g.muK = -decay .* real(sum((Vi.' * Qy) .* V, 2));
end
g.eta(~isfinite(sys.eta)) = 0;
end
