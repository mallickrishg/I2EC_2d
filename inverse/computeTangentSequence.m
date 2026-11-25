function dU = computeTangentSequence(fwd, dTheta, dPhi, varargin)
% COMPUTETANGENTSEQUENCE  Tangent (forward sensitivity) of the integrated states to log viscosity.
%
%   dU = computeTangentSequence(fwd, dTheta)
%   dU = computeTangentSequence(fwd, dTheta, dPhi)
%   dU = computeTangentSequence(fwd, dTheta, dPhi, 'Method', 'grouped')
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
%   'Method'  'separable' (default) or 'grouped', as in
%             computeAdjointGradient: the cost of 'separable' does not
%             depend on the time sampling
%
% OUTPUTS
%   dU      [n x nObs x r] derivative of fwd.U along each direction
%
% NOTES
%   - Tangent equations: d(dy)/dt = A*dy - dTheta.*ydot - dPhi.*decay.*y,
%     jumps -dTheta.*b at every event (observed and characteristic), and
%     the periodicity condition, solved in modal form like the forward.
%   - The sources depend only on the forward states, so they are formed
%     for all segments at once: with the separable form of the pairwise
%     integrals (computePairIntegralSums),
%       sum_j B_ij D_k(lam_i, lam_j, tau_s) y_js
%         = [(B.*C)*(F_k.*Y)](i,s) - F_k(i,s)*[(B.*C)*Y](i,s),
%     plus the close eigenvalue pairs summed directly. The sweep through
%     the segments is then O(n) per segment.
%   - Cost per direction: B = Vi*diag(dTheta)*V (two n x n products) plus
%     a few n x n x (segments) products, whatever the time sampling.
%   - Directions are processed 12 at a time, so memory stays bounded for
%     hundreds of directions (e.g. a full Jacobian for a resolution
%     analysis); the divided differences are recomputed per chunk.
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
% directions in chunks: each holds one or two n x n coupling matrices
% (35 MB each for n = 2103), so hundreds at once would not fit in memory
chunk = 12;
if r > chunk
    dU = zeros(n, numel(fwd.seq.tObs), r);
    for k0 = 1:chunk:r
        kk = k0:min(k0 + chunk - 1, r);
        dU(:, :, kk) = computeTangentSequence(fwd, dTheta(:, kk), dPhi(:, kk), varargin{:});
    end
    return
end
ip = inputParser;
ip.addParameter('Method', 'separable', @(m) any(strcmpi(m, {'separable', 'grouped'})));
ip.parse(varargin{:});
separable = strcmpi(ip.Results.Method, 'separable');
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
if separable
    % all segments at once: separable pairs in two products, close pairs
    % (computeExpDividedDifferences) through a sparse sum
    [C, I, J] = separablePairs(lam, 1e-4);
    F = {shiftedModalIntegrals(lam, segTau, 1), shiftedModalIntegrals(lam, segTau, 2), ...
         shiftedModalIntegrals(lam, segTau, 3)};
    nc = numel(I);
    Psel = sparse(I, (1:nc)', 1, n, nc);
    idx = sub2ind([n n], I, J);
    Dc = {zeros(nc, nall) + 0*lam(1), zeros(nc, nall) + 0*lam(1), zeros(nc, nall) + 0*lam(1)};
    for s = 1:nall
        if segTau(s) == 0, continue, end
        [Dc{1}(:, s), Dc{2}(:, s), Dc{3}(:, s)] = computeExpDividedDifferences(lam(I), lam(J), segTau(s));
    end
    W = lam .* segZ + ch;
    Ch = repmat(ch, 1, nall);
    for k = 1:r
        srcEnd(:, :, k) = -applyPairs(B{k}, W, 1);
        srcInt(:, :, k) = -applyPairs(B{k}, W, 2);
        if doMu
            srcEnd(:, :, k) = srcEnd(:, :, k) - applyPairs(Bm{k}, segZ, 1) - applyPairs(Bm{k}, Ch, 2);
            srcInt(:, :, k) = srcInt(:, :, k) - applyPairs(Bm{k}, segZ, 2) - applyPairs(Bm{k}, Ch, 3);
        end
    end
else
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

    function src = applyPairs(Mx, Y, k)
        % sum_j Mx(i,j)*D_k(lam_i, lam_j, tau_s)*Y(j,s) for all i, s
        MC = Mx .* C;
        src = MC * (F{k} .* Y) - F{k} .* (MC * Y) + Psel * (Mx(idx) .* Dc{k} .* Y(J, :));
    end
end
