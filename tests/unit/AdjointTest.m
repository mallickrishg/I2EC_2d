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

        function verifyGradient(tc, kind, parameters, seq)
            % adjoint gradient against central differences, per parameter
            % (all of them: the small model has a few dozen)
            if nargin < 4, seq = []; end
            tc.setKelvin(kind);
            p0 = tc.parameters();
            fwd = tc.forwardAt(p0, 'auto', seq);
            obs = tc.observations(fwd.sys, fwd);
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
