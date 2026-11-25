function dU = computeTangentSequence(fwd, dTheta, dPhi)
% COMPUTETANGENTSEQUENCE  Tangent (forward sensitivity) of the integrated states to log viscosity.
%
%   dU = computeTangentSequence(fwd, dTheta)
%   dU = computeTangentSequence(fwd, dTheta, dPhi)
%
%   Directional derivatives of fwd.U (solveLinearSequence) for r
%   directions at once, reusing the forward's eigendecomposition. Column k
%   of dTheta perturbs log(sys.eta) of every state row, column k of dPhi
%   log(muK) of every Kelvin row. With a station operator H, H*dU(:,:,k)
%   is the Gauss-Newton Jacobian column for direction k; with r low-rank
%   (e.g. KL) modes the whole Jacobian costs r of these. Build row
%   directions from physical ones with assembleRowPerturbation.
%
% INPUTS
%   fwd     struct from solveLinearSequence
%   dTheta  [n x r] perturbations of log(eta) per state row
%   dPhi    [n x r] perturbations of log(muK) per state row (only Kelvin
%           rows matter)  [zeros]
%
% OUTPUTS
%   dU      [n x nObs x r] derivative of fwd.U along each direction
%
% NOTES
%   - Tangent equations: d(dy)/dt = A*dy - dTheta.*ydot - dPhi.*decay.*y,
%     jumps -dTheta.*b at every event (observed and characteristic), and
%     the periodicity condition, solved in modal form like the forward.
%   - The sources depend only on the forward states, so for each distinct
%     segment length they are formed with one divided-difference matrix
%     (computeExpDividedDifferences) and one product per direction; the
%     sweep through the segments is then O(n) per segment.
%   - Cost per direction: B = Vi*diag(dTheta)*V (two n x n products) plus
%     O(n^2) per distinct segment length.
%   - Adjoint of computeAdjointGradient: for any S,
%     sum(S.*dU) = g.eta'*dTheta + g.muK'*dPhi (tests/unit/AdjointTest).
%
% SEE ALSO
%   solveLinearSequence, computeAdjointGradient, assembleRowPerturbation,
%   computeDisplacementMisfit
%
% AUTHORS
%   Rishav Mallick (EOS), 2026

V = fwd.V; Vi = fwd.Vi; lam = fwd.lam; ch = fwd.ch;
n = numel(lam);
r = size(dTheta, 2);
if nargin < 3 || isempty(dPhi)
    dPhi = zeros(n, r);
end
sys = fwd.sys; seq = fwd.seq;
decay = sys.decay;
dTheta(~isfinite(sys.eta), :) = 0;
doMu = any(decay ~= 0) && any(dPhi(:) ~= 0);
bp = fwd.bp; nbp = numel(bp); nseg = nbp - 1;
tau = diff(bp);
T = seq.T; a = seq.age;

% modal coupling matrices per direction
B = cell(1, r); Bm = cell(1, r);
for k = 1:r
    B{k} = Vi * (dTheta(:, k) .* V);
    if doMu
        Bm{k} = Vi * ((dPhi(:, k) .* decay) .* V);
    end
end

% ---- sources: end value and integral over each segment ----
% segments: [observed..., spin-up to age a, whole cycle]
segTau = [tau, a, T];
segZ = [fwd.zs(:, 1:nseg), fwd.z0, fwd.z0];
nall = numel(segTau);
srcEnd = zeros(n, nall, r) + 0*lam(1);
srcInt = srcEnd;
[ut, grp] = groupLengths(segTau);
for u = 1:numel(ut)
    kk = find(grp == u);
    if ut(u) == 0, continue, end
    [E1, E2, E3] = computeExpDividedDifferences(lam, lam.', ut(u));
    W = lam .* segZ(:, kk) + ch;
    for k = 1:r
        srcEnd(:, kk, k) = -(B{k} .* E1) * W;
        srcInt(:, kk, k) = -(B{k} .* E2) * W;
        if doMu
            srcEnd(:, kk, k) = srcEnd(:, kk, k) - (Bm{k} .* E1) * segZ(:, kk) - (Bm{k} .* E2) * ch;
            srcInt(:, kk, k) = srcInt(:, kk, k) - (Bm{k} .* E2) * segZ(:, kk) - (Bm{k} .* E3) * ch;
        end
    end
end

% ---- sweep per direction ----
[ea] = computeModalIntegrals(lam, a);
dU = zeros(n, numel(seq.tObs), r);
EG = cell(1, nseg);
for i = 1:nseg
    [EG{i}.e, EG{i}.g] = computeModalIntegrals(lam, tau(i));
end
for k = 1:r
    dz0 = (srcEnd(:, nall, k) - B{k} * fwd.b0hat) ./ (-expm1(lam * T));
    dz = ea .* dz0 + srcEnd(:, nall-1, k);
    Bb = B{k} * fwd.bhat;
    duz = zeros(n, 1) + 0*dz;
    DUZ = zeros(n, nbp) + 0*dz;
    for i = 1:nbp
        if i > 1
            duz = duz + EG{i-1}.g .* dz + srcInt(:, i-1, k);
            dz = EG{i-1}.e .* dz + srcEnd(:, i-1, k);
        end
        dz = dz - Bb(:, i);
        DUZ(:, i) = duz;
    end
    dU(:, :, k) = real(V * (DUZ(:, fwd.iObs) - DUZ(:, fwd.iRef)));
end
end
