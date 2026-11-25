classdef (SharedTestFixtures = {I2ECPathFixture}) AdjointTest < matlab.unittest.TestCase
    %ADJOINTTEST  Adjoint and tangent of the linear periodic-plus-sequence
    % forward model (inverse/): gradients of a surface-displacement misfit
    % with respect to log viscosity, Maxwell and Burgers.
    % See docs/dev/adjoint.md.

    properties
        M       % small model
        Evl     % stress kernels (cell-averaged, certified)
        Devl    % displacement kernels at the 7 surface points
        yr = 3.15e7
        Vpl = 1e-9
    end

    methods (TestClassSetup)
        function build(tc)
            tc.M = makeSmallModel(tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder);
            [~, tc.Evl] = evalc('computeAllStressKernelsBem(tc.M.rcv, tc.M.shz, tc.M.bnd)');
            [~, tc.Devl] = evalc('computeAllDisplacementKernelsBem(tc.M.obs, tc.M.rcv, tc.M.shz, tc.M.bnd, 1)');
        end
    end

    methods (TestMethodSetup)
        function setRheology(tc)
            % as TransientRheologyTest, but with spatially varying
            % viscosities so that no two cells are alike
            rcv = tc.M.rcv; shz = tc.M.shz;
            rcv.Vpl = tc.Vpl * ones(rcv.N, 1);
            rcv.faultLaw = 'linear';
            rcv.Asigma = 1e-6*(5e17/sum(rcv.W)) * linspace(0.5, 2, rcv.N)';
            rcv.pinnedPosition = (1:rcv.N)' <= rcv.N/3;
            shz.n = ones(shz.N, 1);
            shz.alpha = 1./(1e-6*1e18*10.^linspace(-0.5, 0.5, shz.N)');
            shz.e22pl = 1e-15 * linspace(-1, 1, shz.N)';
            shz.e23pl = 2e-15 * ones(shz.N, 1);
            shz.etaK = zeros(shz.N, 0);
            shz.muK = zeros(shz.N, 0);
        end
    end

    %% building blocks
    methods (Test)
        function dividedDifferencesMatchMatrixExponential(tc)
            % Opitz: expm(tau*J) of the bidiagonal matrix with the nodes on
            % its diagonal holds the divided differences of exp(x*tau)
            tau = 2.5;
            lam = [-1e3 -3 -1 -0.4 -0.4000001 -0.4+1e-12 -1e-3 -1e-9 0 -1e-14 ...
                   -2+0.5i -2-0.5i 1e-8 -0.39 -50];
            [D1, D2, D3] = computeExpDividedDifferences(lam.', lam, tau);
            err = zeros(1, 3);
            for i = 1:numel(lam)
                for j = 1:numel(lam)
                    a = lam(i); b = lam(j);
                    X = expm(tau*[a 1 0 0; 0 0 1 0; 0 0 b 1; 0 0 0 0]);   % nodes a, 0, b, 0
                    Y = expm(tau*[a 1; 0 b]);
                    ref = [Y(1,2), X(1,3), X(1,4)];
                    err = max(err, abs([D1(i,j), D2(i,j), D3(i,j)] - ref) ./ abs(ref));
                end
            end
            tc.verifyLessThan(err, 1e-13);
        end

        function pairIntegralSumsMatchDirectSum(tc)
            % separable sums over segments against the per-segment sum of
            % computeExpDividedDifferences: eigenvalues over 7 decades, a
            % cluster at every separation (1e-12 to 1e-2 relative), exact
            % repeats, complex pairs; log-spaced segments, a zero length
            % and a near-repeated length
            rng(3);
            base = -10.^(-12 + 7*rand(40, 1));
            clu = -3e-9*(1 + [0; 1e-12; 1e-9; 1e-7; 1e-5; 1e-4; 2e-4; 1e-3; 1e-2]);
            lam = [base; clu; base(1:3); -1e-9 + 2e-10i; -1e-9 - 2e-10i; -5e-7*(1 + 1e-6); -5e-7];
            n = numel(lam);
            tau = [logspace(3, 10, 30), 0, 1e3, 1e3*(1 + 1e-14)];
            ns = numel(tau);
            terms = struct('k', {1, 2, 3}, 'X', {randn(n, ns) + 1i*randn(n, ns), randn(n, ns), randn(n, ns)}, ...
                'Y', {randn(n, ns), randn(n, ns) + 1i*randn(n, ns), randn(n, ns)});
            S = computePairIntegralSums(lam, tau, terms);
            Sref = zeros(n); Sabs = zeros(n);
            for s = 1:ns
                D = cell(1, 3);
                [D{1}, D{2}, D{3}] = computeExpDividedDifferences(lam, lam.', tau(s));
                for t = 1:3
                    P = (terms(t).X(:, s) * terms(t).Y(:, s).') .* D{terms(t).k};
                    Sref = Sref + P; Sabs = Sabs + abs(P);
                end
            end
            tc.verifyLessThan(max(abs(S - Sref) ./ max(Sabs, realmin), [], 'all'), 1e-10);
        end

        function separableMatchesGroupedMethod(tc)
            % adjoint and tangent: the separable sums against one
            % divided-difference matrix per segment length, with
            % log-spaced epochs (every segment length different)
            for kind = ["maxwell", "burgers"]
                tc.setKelvin(kind);
                sys = tc.system();
                seq = tc.sequence(sys);
                seq.tObs = logspace(-2, log10(15), 25)*tc.yr; seq.tRef = seq.tObs(1);
                fwd = solveLinearSequence(sys, seq);
                rng(13);
                n = numel(sys.eta);
                S = randn(size(fwd.U));
                gs = computeAdjointGradient(fwd, S);
                gg = computeAdjointGradient(fwd, S, 'Method', 'grouped');
                tc.verifyEqual(gs.eta, gg.eta, 'AbsTol', 1e-10*max(abs(gg.eta)), char(kind));
                tc.verifyEqual(gs.muK, gg.muK, 'AbsTol', 1e-10*max([abs(gg.muK); realmin]), char(kind));
                dTheta = randn(n, 2); dPhi = randn(n, 2) .* (sys.decay ~= 0);
                dUs = computeTangentSequence(fwd, dTheta, dPhi);
                dUg = computeTangentSequence(fwd, dTheta, dPhi, 'Method', 'grouped');
                tc.verifyEqual(dUs, dUg, 'AbsTol', 1e-10*max(abs(dUg), [], 'all'), char(kind));
            end
        end

        function logSpacedGradientMatchesFiniteDifferences(tc)
            % Burgers, 25 log-spaced epochs from 0.01 yr
            tc.setKelvin("burgers");
            sys = tc.system();
            seq = tc.sequence(sys);
            seq.tObs = logspace(-2, log10(15), 25)*tc.yr; seq.tRef = seq.tObs(1);
            tc.verifyGradient("burgers", 'default', seq);
        end

        function logSpacedGradientConvergesToDenseUniform(tc)
            % a misfit that is an integral over time, J = 1/2*int (h'*U)^2 dt,
            % by the trapezoid rule: on log-spaced epochs (every segment a
            % different length) the gradient converges, at second order,
            % to the one on 4000 uniform epochs
            tc.setKelvin("burgers");
            sys = tc.system();
            rng(21);
            h = randn(1, numel(sys.eta)) * sys.P' * sys.P;   % any fixed functional of the state
            t0 = 0.01*tc.yr; t1 = 15*tc.yr;
            gRef = tc.quadratureGradient(sys, h, linspace(t0, t1, 4000));
            nList = [50 100 200 400 800];
            err = zeros(size(nList));
            for i = 1:numel(nList)
                g = tc.quadratureGradient(sys, h, logspace(log10(t0), log10(t1), nList(i)));
                err(i) = norm(g - gRef)/norm(gRef);
            end
            msg = sprintf('relative error for %s log-spaced epochs: %s', mat2str(nList), mat2str(err, 2));
            tc.verifyLessThan(diff(err), 0, msg);
            tc.verifyLessThan(err(end), 2e-4, msg);
            rate = -diff(log(err))/log(2);
            tc.verifyGreaterThan(rate(end), 1.5, [msg ' (second order expected)']);
        end

        function messySamplingMatchesGroupedMethod(tc)
            % GNSS-like epochs: log-spaced, a daily burst, two long gaps,
            % near-duplicate epochs (1 s, 1 ms and 1 us apart): separable
            % adjoint and tangent against the per-length method
            for kind = ["maxwell", "burgers"]
                tc.setKelvin(kind);
                sys = tc.system();
                seq = tc.sequence(sys);
                seq.tObs = tc.messyEpochs(); seq.tRef = seq.tObs(1);
                fwd = solveLinearSequence(sys, seq);
                rng(17);
                n = numel(sys.eta);
                S = randn(size(fwd.U));
                gs = computeAdjointGradient(fwd, S);
                gg = computeAdjointGradient(fwd, S, 'Method', 'grouped');
                tc.verifyTrue(all(isfinite([gs.eta; gs.muK])), char(kind));
                tc.verifyEqual(gs.eta, gg.eta, 'AbsTol', 1e-10*max(abs(gg.eta)), char(kind));
                tc.verifyEqual(gs.muK, gg.muK, 'AbsTol', 1e-10*max([abs(gg.muK); realmin]), char(kind));
                dTheta = randn(n, 2); dPhi = randn(n, 2) .* (sys.decay ~= 0);
                dUs = computeTangentSequence(fwd, dTheta, dPhi);
                dUg = computeTangentSequence(fwd, dTheta, dPhi, 'Method', 'grouped');
                tc.verifyEqual(dUs, dUg, 'AbsTol', 1e-10*max(abs(dUg), [], 'all'), char(kind));
                % and the tangent stays the adjoint of the gradient
                tc.verifyEqual(sum(S .* dUs(:, :, 1), 'all'), gs.eta'*dTheta(:, 1) + gs.muK'*dPhi(:, 1), ...
                    'RelTol', 1e-10, char(kind));
            end
        end

        function messySamplingAndMissingStationsMatchFiniteDifferences(tc)
            % messy epochs plus stations that start late, stop early, have
            % no data at all, and random dropouts
            tc.setKelvin("burgers");
            sys = tc.system();
            seq = tc.sequence(sys);
            seq.tObs = tc.messyEpochs(); seq.tRef = seq.tObs(1);
            tc.verifyGradient("burgers", 'default', seq, tc.missingPattern(seq.tObs));
        end

        function missingDataCarryNoMisfit(tc)
            % the misfit with NaN entries equals the misfit of the observed
            % entries alone (offsets refitted per row); a row with no data
            % gives no NaN, no misfit and zero Jacobian rows; the adjoint
            % gradient equals Jw'*rw from the tangents
            tc.setKelvin("burgers");
            sys = tc.system();
            seq = tc.sequence(sys);
            seq.tObs = tc.messyEpochs(); seq.tRef = seq.tObs(1);
            fwd = solveLinearSequence(sys, seq);
            obs = tc.observations(sys, fwd);
            miss = tc.missingPattern(seq.tObs);
            obs.d(miss) = NaN;
            n = numel(sys.eta);
            rng(19);
            dTheta = randn(n, 3); dPhi = randn(n, 3) .* (sys.decay ~= 0);
            [J, res] = computeDisplacementMisfit(obs, fwd, computeTangentSequence(fwd, dTheta, dPhi));
            u = obs.H * fwd.U; Jref = 0;
            for i = 1:size(u, 1)
                k = ~miss(i, :);
                if ~any(k), continue, end
                ri = u(i, k) - obs.d(i, k);
                ri = ri - mean(ri);                         % optimal offset (equal weights)
                Jref = Jref + 0.5*obs.W*sum(ri.^2);
            end
            tc.verifyEqual(J, Jref, 'RelTol', 1e-12);
            tc.verifyTrue(all(isfinite(res.S(:))) && all(isfinite(res.Jw(:))));
            dead = all(miss, 2);
            tc.assertTrue(any(dead), 'the pattern has a station with no data');
            Jw = reshape(res.Jw, size(u, 1), [], 3);
            tc.verifyEqual(Jw(dead, :, :), zeros(nnz(dead), size(u, 2), 3));
            tc.verifyEqual(res.offset(dead), zeros(nnz(dead), 1));
            g = computeAdjointGradient(fwd, res.S);
            tc.verifyEqual(res.Jw'*res.rw, dTheta'*g.eta + dPhi'*g.muK, 'RelTol', 1e-10);
        end

        function periodicStateMatchesDirectSolve(tc)
            % modal periodic state of solveLinearSequence against
            % computePeriodicInitialCondition, and the observed window
            % against solveAnalyticIvp chained by hand
            for kind = ["maxwell", "burgers"]
                tc.setKelvin(kind);
                sys = tc.system();
                seq = tc.sequence(sys);
                fwd = solveLinearSequence(sys, seq);
                y0 = computePeriodicInitialCondition(sys.A, sys.c, seq.b0, seq.T);
                tc.verifyEqual(fwd.y0, y0, 'AbsTol', 1e-9*max(abs(y0)), char(kind));

                % chain: start state (age T), event at 0, event at t2
                ys = fwd.ystart + seq.bEvents(:,1);
                t2 = seq.tEvents(2);
                tq = [seq.tRef, seq.tObs];
                U = zeros(size(sys.A, 1), numel(tq));
                [y2, U2] = solveAnalyticIvp(sys.A, sys.c, ys, t2);
                pre = tq <= t2;
                [~, U(:, pre)] = solveAnalyticIvp(sys.A, sys.c, ys, tq(pre));
                [~, Ua] = solveAnalyticIvp(sys.A, sys.c, y2 + seq.bEvents(:,2), tq(~pre) - t2);
                U(:, ~pre) = U2 + Ua;
                Uref = U(:, 2:end) - U(:, 1);
                tc.verifyEqual(fwd.U, Uref, 'AbsTol', 1e-8*max(abs(Uref), [], 'all'), char(kind));
            end
        end
    end

    %% gradients against finite differences
    methods (Test)
        function maxwellGradientMatchesFiniteDifferences(tc)
            tc.verifyGradient("maxwell", 'default');
        end

        function burgersGradientMatchesFiniteDifferences(tc)
            tc.verifyGradient("burgers", 'default');
        end

        function burgersRawGradientMatchesFiniteDifferences(tc)
            tc.verifyGradient("burgers", 'raw');
        end

        function periodicityAdjointAlone(tc)
            % a misfit that sees only the spin-up: observations before any
            % observed event, sources only from the start state onwards
            sys = tc.system();
            seq = tc.sequence(sys);
            seq.tEvents = []; seq.bEvents = [];
            seq.age = 0.4*seq.T;
            tc.verifyGradient("maxwell", 'default', seq);
        end
    end

    %% tangent
    methods (Test)
        function tangentPassesDotProductTest(tc)
            for kind = ["maxwell", "burgers"]
                tc.setKelvin(kind);
                sys = tc.system();
                fwd = solveLinearSequence(sys, tc.sequence(sys));
                rng(7);
                n = numel(sys.eta);
                dTheta = randn(n, 3); dPhi = randn(n, 3) .* (sys.decay ~= 0);
                S = randn(size(fwd.U));
                dU = computeTangentSequence(fwd, dTheta, dPhi);
                g = computeAdjointGradient(fwd, S);
                for k = 1:3
                    lhs = sum(S .* dU(:, :, k), 'all');
                    rhs = g.eta' * dTheta(:, k) + g.muK' * dPhi(:, k);
                    tc.verifyEqual(lhs, rhs, 'RelTol', 1e-10, sprintf('%s, direction %d', kind, k));
                end
            end
        end

        function tangentMatchesFiniteDifferences(tc)
            tc.setKelvin("burgers");
            sys = tc.system();
            p0 = tc.parameters();
            rng(11);
            dp = tc.randomDirection(p0);
            fwd = solveLinearSequence(sys, tc.sequence(sys));
            [dTheta, dPhi] = assembleRowPerturbation(sys, dp);
            dU = computeTangentSequence(fwd, dTheta, dPhi);
            h = 1e-4;
            Up = tc.forwardAt(tc.shift(p0, dp, h), fwd.balance);
            Um = tc.forwardAt(tc.shift(p0, dp, -h), fwd.balance);
            tc.setParameters(p0);
            dUfd = (Up.U - Um.U) / (2*h);
            tc.verifyEqual(dU, dUfd, 'AbsTol', 1e-6*max(abs(dUfd), [], 'all'));
        end

        function parameterMapsAreTransposes(tc)
            tc.setKelvin("burgers");
            sys = tc.system();
            rng(5);
            n = numel(sys.eta);
            g.eta = randn(n, 1); g.muK = randn(n, 1) .* (sys.decay ~= 0);
            for kind = {'default', 'raw'}
                gp = extractViscosityGradient(g, sys, 'Parameters', kind{1});
                dp = structfun(@(x) randn(size(x)), gp, 'UniformOutput', false);
                [dTheta, dPhi] = assembleRowPerturbation(sys, dp, 'Parameters', kind{1});
                lhs = sum(cellfun(@(f) sum(gp.(f) .* dp.(f), 'all'), fieldnames(gp)));
                tc.verifyEqual(lhs, g.eta'*dTheta + g.muK'*dPhi, 'RelTol', 1e-12, kind{1});
            end
        end

        function tangentInChunksMatchesOneAtATime(tc)
            % computeTangentSequence processes directions 12 at a time
            tc.setKelvin("burgers");
            sys = tc.system();
            fwd = solveLinearSequence(sys, tc.sequence(sys));
            rng(9);
            n = numel(sys.eta);
            dTheta = randn(n, 15); dPhi = randn(n, 15) .* (sys.decay ~= 0);
            dU = computeTangentSequence(fwd, dTheta, dPhi);
            for k = [1 12 13 15]
                dUk = computeTangentSequence(fwd, dTheta(:, k), dPhi(:, k));
                tc.verifyEqual(dU(:, :, k), dUk, 'AbsTol', 1e-12*max(abs(dUk), [], 'all'), sprintf('direction %d', k));
            end
        end

        function balancedModesAreAccurateForSpreadViscosity(tc)
            % viscosities spread over several decades: the balanced modes
            % (Vi from the balanced eigenvectors, not from inverting V)
            % match the matrix exponential; both variants agree with each
            % other on the default model
            sys = tc.system();
            seq = tc.sequence(sys);
            fU = solveLinearSequence(sys, seq, 'Balance', false);
            fB = solveLinearSequence(sys, seq, 'Balance', true);
            tc.verifyEqual(fB.U, fU.U, 'AbsTol', 1e-8*max(abs(fU.U), [], 'all'));

            rng(4);
            shz = tc.M.shz;
            shz.alpha = 1./(1e-6*10.^(19 + 1.5*randn(shz.N, 1)));
            shz.etaK = 1e-6*10.^(18 + 1.5*randn(shz.N, 1)); shz.muK = 1.5e4*ones(shz.N, 1);
            sys = tc.system();
            seq = tc.sequence(sys);
            fB = solveLinearSequence(sys, seq, 'Balance', true);
            n = size(sys.A, 1);
            X = expm([sys.A, sys.c; zeros(1, n+1)] * seq.T);
            yT = X(1:n, 1:n)*fB.y0 + X(1:n, end);            % one cycle from the periodic state
            tc.verifyEqual(yT + seq.b0, fB.y0, 'AbsTol', 1e-6*max(abs(fB.y0)));
        end

        function rowsWithoutViscosityGetNoGradient(tc)
            % gravity readiness: a row with eta = NaN (no viscosity) is
            % carried by the solve but is not a parameter
            sys = tc.system();
            sys.eta(1) = NaN;
            fwd = solveLinearSequence(sys, tc.sequence(sys));
            g = computeAdjointGradient(fwd, randn(size(fwd.U)));
            tc.verifyEqual(g.eta(1), 0);
            tc.verifyTrue(all(isfinite(g.eta)));
        end
    end

    methods (Access = private)
        function setKelvin(tc, kind)
            shz = tc.M.shz;
            switch kind
                case "maxwell"
                    shz.etaK = zeros(shz.N, 0); shz.muK = zeros(shz.N, 0);
                case "burgers"     % tauK ~ 0.3 - 3 yr, varying per cell
                    shz.etaK = 1e-6 * 1e18 * 10.^linspace(0.5, -0.5, shz.N)';
                    shz.muK = 3e4 * linspace(0.8, 1.2, shz.N)';
            end
        end

        function sys = system(tc)
            sys = assembleLinearSystem(tc.Evl, tc.M.rcv, tc.M.shz, tc.M.rcv.pinnedPosition, ...
                tc.events(tc.M.rcv.pinnedPosition));
        end

        function sc = events(tc, locked)
            % event 1: 200 yr of slip deficit on the locked patches;
            % event 2: half of it on the shallower half
            evl = tc.Evl;
            slip = [200*tc.yr*tc.Vpl*locked, 100*tc.yr*tc.Vpl*(locked & (1:numel(locked))' <= 2)];
            sc = struct('Nevents', 2, 'Timing', [0 0], 'dtau', evl.KK*slip, ...
                'dsigma22', evl.KL(:,:,1)*slip, 'dsigma23', evl.KL(:,:,2)*slip);
        end

        function seq = sequence(tc, sys)
            % spin-up: event 1 every 200 yr; observed: event 1 at t = 0,
            % event 2 at 6 yr; 14 epochs from 0.3 to 15 yr, reference at 0.1 yr
            yr = tc.yr;
            seq = struct('T', 200*yr, 'b0', sys.b(:,1), 'age', 200*yr, ...
                'tEvents', [0 6*yr], 'bEvents', sys.b(:, [1 2]), ...
                'tObs', linspace(0.3, 15, 14)*yr, 'tRef', 0.1*yr);
        end

        function t = messyEpochs(tc)
            % log-spaced epochs over 15 yr with a daily burst for 30 days
            % at 0.5 yr, gaps over 2-6 yr and 9-11 yr, near-duplicates 1 s,
            % 1 ms and 1 us after three epochs, and a last lone epoch
            yr = tc.yr;
            t = logspace(-2, log10(12), 60)*yr;
            t = t(~(t > 2*yr & t < 6*yr) & ~(t > 9*yr & t < 11*yr));
            t = [t, 0.5*yr + (0:29)*yr/365, t(10) + 1, t(20) + 1e-3, t(30) + 1e-6, 15*yr];
            t = unique(t);
        end

        function miss = missingPattern(tc, t)
            % 7 stations x 2 components: stations 1-2 start at 5 yr,
            % station 3 stops at 3 yr, station 4 has no data, 10% random
            % dropouts elsewhere
            rng(23);
            ns = size(tc.M.obs, 1);
            m1 = rand(ns, numel(t)) < 0.1;
            m1(1:2, :) = m1(1:2, :) | t < 5*tc.yr;
            m1(3, :) = m1(3, :) | t > 3*tc.yr;
            m1(4, :) = true;
            miss = [m1; m1];                                 % both components
        end

        function g = quadratureGradient(tc, sys, h, t)
            % gradient of J = 1/2*sum_j w_j*(h*U_j)^2, trapezoid weights w on t
            seq = tc.sequence(sys);
            seq.tObs = t; seq.tRef = t(1);
            fwd = solveLinearSequence(sys, seq);
            wq = [diff(t), 0]/2 + [0, diff(t)]/2;
            u = h * fwd.U;
            gr = computeAdjointGradient(fwd, h' * (wq .* u));
            g = gr.eta;
        end

        function obs = observations(tc, sys, fwd)
            % both components at the 7 stations; data: the model's
            % prediction plus a smooth perturbation, so residuals are O(1)
            obs.H = assembleStationKernels(tc.Devl, sys);
            u = obs.H * fwd.U;
            s = max(abs(u), [], 'all');
            obs.d = u + 0.1*s*sin((1:size(u, 1))' * (1:size(u, 2)) / 3);
            obs.W = 1/(0.01*s)^2;
            obs.offsets = true;
        end

        function p = parameters(tc)
            rcv = tc.M.rcv; shz = tc.M.shz;
            p.fault = log(rcv.Asigma);
            p.etaM = log(1./shz.alpha);
            [etaK, muK] = extractKelvinElements(shz);
            if ~isempty(etaK)
                p.tauK = log(etaK ./ muK);
                p.muKrel = log(muK / 3e4);
            end
        end

        function setParameters(tc, p)
            tc.M.rcv.Asigma = exp(p.fault);
            tc.M.shz.alpha = exp(-p.etaM);
            if isfield(p, 'tauK')
                tc.M.shz.muK = 3e4 * exp(p.muKrel);
                tc.M.shz.etaK = exp(p.tauK) .* tc.M.shz.muK;
            end
        end

        function dp = randomDirection(~, p)
            dp = structfun(@(x) randn(size(x)), p, 'UniformOutput', false);
        end

        function q = shift(~, p, dp, h)
            q = p;
            for f = fieldnames(p)'
                q.(f{1}) = p.(f{1}) + h*dp.(f{1});
            end
        end

        function fwd = forwardAt(tc, p, balance, seq)
            tc.setParameters(p);
            sys = tc.system();
            if nargin < 4 || isempty(seq)
                seq = tc.sequence(sys);
            else            % same time line, jumps of this system
                seq.b0 = sys.b(:,1);
                if ~isempty(seq.tEvents), seq.bEvents = sys.b(:, 1:numel(seq.tEvents)); end
            end
            fwd = solveLinearSequence(sys, seq, 'Balance', balance);
        end

        function verifyGradient(tc, kind, parameters, seq, missing)
            % adjoint gradient against central differences, per parameter
            % (all of them: the small model has a few dozen)
            if nargin < 4, seq = []; end
            tc.setKelvin(kind);
            p0 = tc.parameters();
            fwd = tc.forwardAt(p0, 'auto', seq);
            obs = tc.observations(fwd.sys, fwd);
            if nargin > 4 && ~isempty(missing)
                obs.d(missing) = NaN;
            end
            [~, res] = computeDisplacementMisfit(obs, fwd);
            g = computeAdjointGradient(fwd, res.S);
            gp = extractViscosityGradient(g, fwd.sys, 'Parameters', parameters);

            if strcmp(parameters, 'raw')            % FD in raw logs
                toRaw = @(p) rawOf(p);
                fromRaw = @(q) defaultOf(q);
            else
                toRaw = @(p) p; fromRaw = @(q) q;
            end
            q0 = toRaw(p0);
            h = 1e-4;
            Jat = @(q) computeDisplacementMisfit(obs, tc.forwardAt(fromRaw(q), fwd.balance, seq));
            fields = fieldnames(gp)';
            gfd = gp;
            for f = fields
                gfd.(f{1}) = zeros(size(gp.(f{1})));
                for i = 1:numel(gfd.(f{1}))
                    if strcmp(f{1}, 'fault') && tc.M.rcv.pinnedPosition(i), continue, end
                    qp = q0; qp.(f{1})(i) = qp.(f{1})(i) + h;
                    qm = q0; qm.(f{1})(i) = qm.(f{1})(i) - h;
                    gfd.(f{1})(i) = (Jat(qp) - Jat(qm)) / (2*h);
                end
            end
            % tolerance relative to the largest component of the whole
            % gradient: weakly sensitive parameters are at the FD noise
            scale = max(cellfun(@(f) max(abs(gfd.(f)(:))), fields));
            for f = fields
                tc.verifyEqual(gp.(f{1}), gfd.(f{1}), 'AbsTol', 1e-6*scale, ...
                    sprintf('%s, %s, %s', kind, parameters, f{1}));
            end
            tc.setParameters(p0);

            function q = rawOf(p)
                q = p;
                if isfield(p, 'tauK')
                    q = rmfield(q, {'tauK', 'muKrel'});
                    q.muK = p.muKrel + log(3e4);
                    q.etaK = p.tauK + q.muK;
                end
            end
            function p = defaultOf(q)
                p = q;
                if isfield(q, 'etaK')
                    p = rmfield(p, {'etaK', 'muK'});
                    p.tauK = q.etaK - q.muK;
                    p.muKrel = q.muK - log(3e4);
                end
            end
        end
    end
end
