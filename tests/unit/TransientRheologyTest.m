classdef (SharedTestFixtures = {I2ECPathFixture}) TransientRheologyTest < matlab.unittest.TestCase
    %TRANSIENTRHEOLOGYTEST  Burgers and Andrade (Kelvin elements in series
    % with the Maxwell dashpot) in both solvers: assembleLinearSystem for
    % the semi-analytical path, odeViscoelastic for the numerical one,
    % computeAndradePronySeries for Andrade. See docs/dev/transient-rheology.md.

    properties
        M       % small model
        Evl     % default kernels (cell-averaged, certified)
        yr = 3.15e7
        Vpl = 1e-9
    end

    methods (TestClassSetup)
        function build(tc)
            tc.M = makeSmallModel(tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder);
            [~, tc.Evl] = evalc('computeAllStressKernelsBem(tc.M.rcv, tc.M.shz, tc.M.bnd)');
        end
    end

    methods (TestMethodSetup)
        function setRheology(tc)
            % linear viscous fault, shallow third locked, 1e18 Pa s mantle,
            % no Kelvin elements (rcv and shz are handles: reset per test)
            rcv = tc.M.rcv; shz = tc.M.shz;
            rcv.Vpl = tc.Vpl * ones(rcv.N, 1);
            rcv.faultLaw = 'linear';
            rcv.Asigma = 1e-6*(5e17/sum(rcv.W)) * ones(rcv.N, 1);
            rcv.pinnedPosition = (1:rcv.N)' <= rcv.N/3;
            shz.n = ones(shz.N, 1);
            shz.alpha = 1/(1e18*1e-6) * ones(shz.N, 1);
            shz.e22pl = 1e-15 * linspace(-1, 1, shz.N)';
            shz.e23pl = 2e-15 * ones(shz.N, 1);
            shz.etaK = zeros(shz.N, 0);
            shz.muK = zeros(shz.N, 0);
        end
    end

    %% Maxwell is the default and is unchanged
    methods (Test)
        function maxwellSystemIsTheExistingOne(tc)
            % no Kelvin elements: A, c, b are the expressions of the
            % Maxwell solver (newtonianCycle) bit for bit
            rcv = tc.M.rcv; shz = tc.M.shz; evl = tc.Evl;
            locked = rcv.pinnedPosition;
            sc = tc.event(locked);
            sys = assembleLinearSystem(evl, rcv, shz, locked, sc);
            eta = [rcv.Asigma(~locked); 1./shz.alpha; 1./shz.alpha];
            A = assembleStressKernel(evl, locked) ./ eta;
            vlock = rcv.Vpl .* locked;
            c = -[evl.KK(~locked,:)*vlock; evl.KL(:,:,1)*vlock; evl.KL(:,:,2)*vlock] ./ eta ...
                - A*[rcv.Vpl(~locked); shz.e22pl; shz.e23pl];
            b = [sc.dtau(~locked); sc.dsigma22; sc.dsigma23] ./ eta;
            tc.verifyEqual(sys.A, A);
            tc.verifyEqual(sys.c, c);
            tc.verifyEqual(sys.b, b);
            tc.verifyEqual(sys.p, 0);
            tc.verifyEqual(full(sys.P), eye(numel(c)));
        end

        function infiniteKelvinViscosityIsMaxwell(tc)
            % etaK = Inf: no Kelvin state, the Maxwell system exactly; the
            % ODE right-hand side has zero Kelvin rate and the Maxwell rest
            rcv = tc.M.rcv; shz = tc.M.shz; evl = tc.Evl;
            locked = rcv.pinnedPosition;
            ref = assembleLinearSystem(evl, rcv, shz, locked, tc.event(locked));
            Y = tc.someOdeState(0);
            Yp0 = odeViscoelastic(0, Y, rcv, shz, evl);
            shz.etaK = Inf(shz.N, 2); shz.muK = 3e4 * ones(shz.N, 2);
            sys = assembleLinearSystem(evl, rcv, shz, locked, tc.event(locked));
            tc.verifyEqual(sys.A, ref.A);
            tc.verifyEqual(sys.c, ref.c);
            tc.verifyEqual(sys.b, ref.b);
            Y = tc.someOdeState(2);
            Yp = odeViscoelastic(0, Y, rcv, shz, evl);
            tc.verifyEqual(Yp(1:numel(Yp0)), Yp0);
            tc.verifyEqual(Yp(numel(Yp0)+1:end), zeros(4*shz.N, 1));
        end

        function largeKelvinViscosityConvergesToMaxwell(tc)
            % Burgers with etaK -> Inf: the periodic cycle tends to the
            % Maxwell cycle, with an error that falls like 1/etaK
            rcv = tc.M.rcv; shz = tc.M.shz; evl = tc.Evl;
            locked = rcv.pinnedPosition; T = 100*tc.yr; tvec = [0.1 1 10 99]*tc.yr;
            ref = tc.periodicRates(evl, rcv, shz, locked, T, tvec);
            etaKs = 1e-6 * [1e20 1e22 1e24];      % MPa s
            err = zeros(size(etaKs));
            for i = 1:numel(etaKs)
                shz.etaK = etaKs(i) * ones(shz.N, 1); shz.muK = 3e4 * ones(shz.N, 1);
                X = tc.periodicRates(evl, rcv, shz, locked, T, tvec);
                err(i) = max(abs(X - ref), [], 'all') / max(abs(ref), [], 'all');
            end
            tc.verifyLessThan(err(end), 1e-5, sprintf('errors %s', mat2str(err, 2)));
            tc.verifyLessThan(err(2:end) ./ err(1:end-1), 0.1, 'error falls like 1/etaK');
        end
    end

    %% Single cell: closed-form Burgers responses
    methods (Test)
        function singleCellBurgersCreep(tc)
            % constant stress (no elastic feedback): strain
            %   e(t) = s0*(t/etaM + (1 - exp(-t/tauK))/muK)
            [rcv, shz, evl, s0, etaM, etaK, muK] = tc.singleCell(0);
            tauK = etaK/muK; t = [0.01 0.3 1 3 10 30]*tc.yr;
            eAn = s0*(t/etaM + (1 - exp(-t/tauK))/muK);
            edotAn = s0*(1/etaM + exp(-t/tauK)/etaK);

            % semi-analytical: rates [eM22; eM23; eK22; eK23], integral = strain
            sys = assembleLinearSystem(evl, rcv, shz, true);
            y0 = [s0/etaM; 0; s0/etaK; 0];
            [Y, U] = solveAnalyticIvp(sys.A, sys.c, y0, t);
            [~, e22dot] = extractRatesFromLinearState(Y, sys);
            [~, e22, e23] = extractRatesFromLinearState(U, sys);
            tc.verifyEqual(e22, eAn, 'RelTol', 1e-10);
            tc.verifyEqual(e23, zeros(size(t)), 'AbsTol', 1e-20);
            tc.verifyEqual(e22dot, edotAn, 'RelTol', 1e-10);

            % numerical: state [fault; s22; s23; sK22; sK23], sK = muK*epsK
            [Ynum, edotNum] = tc.integrateSingleCell(rcv, shz, evl, [0; s0; 0; 0; 0], t);
            tc.verifyEqual(Ynum(:,2)', s0*ones(size(t)), 'RelTol', 1e-12, 'stress stays constant');
            tc.verifyEqual(Ynum(:,4)'/muK, s0/muK*(1 - exp(-t/tauK)), 'RelTol', 1e-7, 'Kelvin strain');
            tc.verifyEqual(edotNum, edotAn, 'RelTol', 1e-7);
        end

        function singleCellBurgersRelaxation(tc)
            % fixed total strain, self-stiffness G: ds/dt = -G*edot, so
            %   s(t) = a1 exp(r1 t) + a2 exp(r2 t),
            %   r^2 + r (G/etaM + (G+muK)/etaK) + G muK/(etaM etaK) = 0,
            %   s(0) = s0, s'(0) = -G s0 (1/etaM + 1/etaK)
            G = 3e4;
            [rcv, shz, evl, s0, etaM, etaK, muK] = tc.singleCell(G);
            bq = G/etaM + (G + muK)/etaK; cq = G*muK/(etaM*etaK);
            r = (-bq + [1 -1]*sqrt(bq^2 - 4*cq))/2;
            ds0 = -G*s0*(1/etaM + 1/etaK);
            a1 = (ds0 - r(2)*s0)/(r(1) - r(2)); a2 = s0 - a1;
            t = [0.01 0.3 1 3 10 30]*tc.yr;
            sAn = a1*exp(r(1)*t) + a2*exp(r(2)*t);
            % Kelvin strain from s = muK epsK + etaK depsK/dt, epsK(0) = 0
            eKAn = (a1/(etaK*r(1) + muK))*(exp(r(1)*t) - exp(-muK/etaK*t)) + ...
                   (a2/(etaK*r(2) + muK))*(exp(r(2)*t) - exp(-muK/etaK*t));

            sys = assembleLinearSystem(evl, rcv, shz, true);
            Y = solveAnalyticIvp(sys.A, sys.c, [s0/etaM; 0; s0/etaK; 0], t);
            tc.verifyEqual(etaM*Y(1,:), sAn, 'RelTol', 1e-9, 'stress, semi-analytical');
            tc.verifyEqual((etaM*Y(1,:) - etaK*Y(3,:))/muK, eKAn, 'RelTol', 1e-8, 'Kelvin strain, semi-analytical');

            Ynum = tc.integrateSingleCell(rcv, shz, evl, [0; s0; 0; 0; 0], t);
            tc.verifyEqual(Ynum(:,2)', sAn, 'RelTol', 1e-7, 'stress, numerical');
            tc.verifyEqual(Ynum(:,4)'/muK, eKAn, 'RelTol', 1e-7, 'Kelvin strain, numerical');
        end

        function andradeCreepFollowsTheOneThirdPowerLaw(tc)
            % Prony series of Kelvin elements creeps like beta*s0*t^(1/3)
            % on [tmin, tmax], within the reported fit error
            [rcv, shz, evl, s0, etaM] = tc.singleCell(0);
            beta = 1e-7; tmin = 1e-3*tc.yr; tmax = 2500*tc.yr;
            fit = computeAndradePronySeries(beta, tmin, tmax);
            tc.verifyLessThan(fit.maxRelErr, 1e-2, 'default fit within 1%');
            fit2 = computeAndradePronySeries(beta, tmin, tmax, 'TermsPerDecade', 2, 'Margin', 1);
            tc.verifyLessThan(fit2.maxRelErr, 1e-4, 'finer fit within 1e-4');
            shz.etaK = fit.etaK; shz.muK = fit.muK;
            sys = assembleLinearSystem(evl, rcv, shz, true);
            y0 = [s0/etaM; 0; zeros(2*sys.p, 1)];
            y0(vertcat(sys.kelvin.rows22)) = s0 ./ fit.etaK';
            t = logspace(log10(tmin), log10(tmax), 30);
            [~, U] = solveAnalyticIvp(sys.A, sys.c, y0, t);
            [~, e22] = extractRatesFromLinearState(U, sys);
            transient = e22 - s0*t/etaM;
            tc.verifyEqual(transient, beta*s0*t.^(1/3), 'RelTol', 1.01*fit.maxRelErr);
        end
    end

    %% Coupled model: spectrum, cycle mean, both solvers agree
    methods (Test)
        function burgersAndAndradeHaveNoGrowingModes(tc)
            % springs and dashpots: real, negative spectrum; with the slip
            % deficit of a cycle released (slip = Vpl*T on the locked zone)
            % the cycle mean of the slip and strain rates is the long-term
            % rate, and the Kelvin rates integrate to zero
            rcv = tc.M.rcv; shz = tc.M.shz; evl = tc.Evl;
            locked = rcv.pinnedPosition; T = 100*tc.yr;
            slip = tc.Vpl*T*locked;
            sc = struct('Nevents', 1, 'Timing', 0, 'dtau', evl.KK*slip, ...
                'dsigma22', evl.KL(:,:,1)*slip, 'dsigma23', evl.KL(:,:,2)*slip);
            for kind = ["burgers" "andrade"]
                tc.setKelvin(kind, T);
                sys = assembleLinearSystem(evl, rcv, shz, locked, sc);
                lam = eig(sys.A);
                tc.verifyLessThan(max(real(lam)), 0, kind + ": no growing modes");
                tc.verifyLessThan(max(abs(imag(lam))), 1e-6*max(abs(lam)), kind + ": real spectrum");
                y0 = computePeriodicInitialCondition(sys.A, sys.c, sys.b, T);
                [~, U] = solveAnalyticIvp(sys.A, sys.c, y0, T);
                [Vm, e22m, e23m, parts] = extractRatesFromLinearState(U/T, sys);
                tc.verifyEqual(Vm(~locked), rcv.Vpl(~locked), 'RelTol', 1e-7, kind + ": mean slip rate");
                emax = max(abs([shz.e22pl; shz.e23pl]));
                tc.verifyEqual([e22m; e23m], [shz.e22pl; shz.e23pl], 'AbsTol', 1e-7*emax, kind + ": mean strain rate");
                tc.verifyEqual([parts.e22K(:); parts.e23K(:)], zeros(2*numel(parts.e22K), 1), ...
                    'AbsTol', 1e-7*emax, kind + ": Kelvin rates integrate to zero");
            end
        end

        function burgersPeriodicCycleMatchesSpunUpNumericalCycle(tc)
            tc.setKelvin("burgers", 100*tc.yr);
            tc.verifySpinUpMatchesPeriodicCycle(100*tc.yr);
        end

        function andradePeriodicCycleMatchesSpunUpNumericalCycle(tc)
            tc.setKelvin("andrade", 100*tc.yr);
            tc.verifySpinUpMatchesPeriodicCycle(100*tc.yr);
        end

        function powerLawBurgersStartsInSteadyState(tc)
            % numerical solver only: n = 3 Maxwell element with Kelvin
            % elements; the long-term state (relaxed Kelvin elements) does
            % not change without earthquakes. Stress errors at the ode45
            % tolerance (1e-10 MPa) give strain-rate errors of 1e-10/etaK
            % = 1e-21 1/s, 1e-6 of the rates here.
            rcv = tc.M.rcv; shz = tc.M.shz;
            shz.n = 3 * ones(shz.N, 1);
            rcv.faultLaw = 'ratestrengthening'; rcv.Asigma = 0.5 * ones(rcv.N, 1);
            rcv.pinnedPosition = false(rcv.N, 1);
            shz.etaK = 1e-6 * [1e17 1e18] .* ones(shz.N, 2); shz.muK = 3e4 * ones(shz.N, 2);
            sc = struct('Nevents', 1, 'Timing', 10*tc.yr, ...
                'dtau', zeros(rcv.N,1), 'dsigma22', zeros(shz.N,1), 'dsigma23', zeros(shz.N,1));
            [~, ~, V, e22, e23, Y] = evalc('runImposedEarthquakeCycles(rcv, shz, tc.Evl, sc, 1, 50*tc.yr)');
            tc.verifySize(Y, [size(V,1), rcv.N + 2*shz.N + 4*shz.N]);
            tc.verifyEqual(V, repmat(rcv.Vpl', size(V,1), 1), 'RelTol', 1e-8);
            tc.verifyEqual(e22, repmat(shz.e22pl', size(V,1), 1), 'AbsTol', 1e-5*1e-15);
            tc.verifyEqual(e23, repmat(shz.e23pl', size(V,1), 1), 'AbsTol', 1e-5*2e-15);
            % the stiff solver option gives the same
            [~, ~, V2, e22b] = evalc('runImposedEarthquakeCycles(rcv, shz, tc.Evl, sc, 1, 50*tc.yr, @ode15s)');
            tc.verifyEqual(V2, repmat(rcv.Vpl', size(V2,1), 1), 'RelTol', 1e-8);
            tc.verifyEqual(e22b, repmat(shz.e22pl', size(V2,1), 1), 'AbsTol', 1e-5*1e-15);
        end

        function invalidKelvinElementsAreRejected(tc)
            shz = tc.M.shz;
            shz.etaK = ones(shz.N, 1); shz.muK = zeros(shz.N, 1);
            tc.verifyError(@() extractKelvinElements(shz), 'extractKelvinElements:invalid');
            shz.etaK = -ones(shz.N, 1); shz.muK = ones(shz.N, 1);
            tc.verifyError(@() extractKelvinElements(shz), 'extractKelvinElements:invalid');
            shz.etaK = ones(shz.N, 2); shz.muK = ones(shz.N, 1);
            tc.verifyError(@() extractKelvinElements(shz), 'extractKelvinElements:invalid');
            shz.etaK = [1 Inf]; shz.muK = [2 0];             % one row for all cells
            [etaK, muK, p] = extractKelvinElements(shz);
            tc.verifyEqual(p, 2);
            tc.verifyEqual(etaK, repmat([1 Inf], shz.N, 1));
            tc.verifyEqual(muK, repmat([2 0], shz.N, 1));
            shz.n = 3 * ones(shz.N, 1);
            tc.verifyError(@() assembleLinearSystem(tc.Evl, tc.M.rcv, shz, false(tc.M.rcv.N, 1)), ...
                'assembleLinearSystem:nonlinear');
        end
    end

    methods (Access = private)
        function setKelvin(tc, kind, T)
            % Burgers: one element, tauK ~ 1 yr; Andrade: Prony series on
            % [1 yr, 0.3 T] (short enough to spin up in a test)
            shz = tc.M.shz;
            switch kind
                case "burgers"
                    shz.etaK = 1e-6 * 1e18 * ones(shz.N, 1);
                    shz.muK = 3e4 * ones(shz.N, 1);
                case "andrade"
                    fit = computeAndradePronySeries(1e-7 * ones(shz.N, 1), tc.yr, 0.3*T);
                    shz.etaK = fit.etaK; shz.muK = fit.muK;
            end
        end

        function sc = event(tc, locked, t1)
            % 2 m of slip on the locked patches
            if nargin < 3, t1 = 0; end
            evl = tc.Evl; slip = 2 * locked;
            sc = struct('Nevents', 1, 'Timing', t1, 'dtau', evl.KK*slip, ...
                'dsigma22', evl.KL(:,:,1)*slip, 'dsigma23', evl.KL(:,:,2)*slip);
        end

        function X = periodicRates(tc, evl, rcv, shz, locked, T, tvec)
            sys = assembleLinearSystem(evl, rcv, shz, locked, tc.event(locked));
            y0 = computePeriodicInitialCondition(sys.A, sys.c, sys.b, T);
            [V, e22, e23] = extractRatesFromLinearState(solveAnalyticIvp(sys.A, sys.c, y0, tvec), sys);
            X = [V(~locked,:); e22; e23];
        end

        function Y = someOdeState(tc, p)
            % arbitrary numerical state with p Kelvin columns
            rcv = tc.M.rcv; shz = tc.M.shz;
            rng(3);
            Y = [tc.Vpl*(1 + rand(rcv.N, 1)); 0.1*randn(2*shz.N, 1); 1e-6*randn(2*shz.N*p, 1)];
        end

        function [rcv, shz, evl, s0, etaM, etaK, muK] = singleCell(~, G)
            % one locked patch (no fault coupling) and one mantle cell with
            % self-stiffness -G in both components; Burgers 1e19 / 1e18 Pa s
            etaM = 1e13; etaK = 1e12; muK = 3e4; s0 = 1;          % MPa s, MPa
            rcv = struct('N', 1, 'dgf', 1, 'Vpl', 0, 'Asigma', 1, 'pinnedPosition', true);
            shz = struct('N', 1, 'dgf', 2, 'alpha', 1/etaM, 'n', 1, ...
                'e22pl', 0, 'e23pl', 0, 'etaK', etaK, 'muK', muK);
            LL = zeros(1, 1, 2, 2); LL(1,1,1,1) = -G; LL(1,1,2,2) = -G;
            evl = struct('KK', 0, 'KL', zeros(1,1,2), 'LK', zeros(1,1,2), 'LL', LL);
        end

        function [Ynum, edot22] = integrateSingleCell(~, rcv, shz, evl, Y0, t)
            opt = odeset('AbsTol', 1e-14, 'RelTol', 1e-10);
            [~, Ynum] = ode45(@(tt, y) odeViscoelastic(tt, y, rcv, shz, evl), [0 t], Y0, opt);
            Ynum = Ynum(2:end, :);
            [~, edot22] = extractRatesFromOdeState(Ynum, rcv, shz);
            edot22 = edot22';
        end

        function verifySpinUpMatchesPeriodicCycle(tc, Trecur)
            % numerical spin-up (runImposedEarthquakeCycles, ode45) against
            % the semi-analytical periodic cycle, slip and strain rates
            rcv = tc.M.rcv; shz = tc.M.shz; evl = tc.Evl;
            locked = rcv.pinnedPosition; t1 = 30*tc.yr;
            sc = tc.event(locked, t1);
            sys = assembleLinearSystem(evl, rcv, shz, locked, sc);
            y0 = computePeriodicInitialCondition(sys.A, sys.c, sys.b, Trecur);

            rho = max(abs(exp(eig(sys.A)*Trecur)));
            nreps = ceil(log(1e-8)/log(rho)) + 1;
            tc.assertLessThan(nreps, 30, sprintf('spin-up too slow: rho = %.3f', rho));
            [~, t, V, e22, e23] = evalc('runImposedEarthquakeCycles(rcv, shz, evl, sc, nreps, Trecur)');

            post = [false; diff(t) == 0];                 % post-event copy of t1
            tau = t - t1 + Trecur*(t < t1 | (t == t1 & ~post));
            [Van, e22an, e23an] = extractRatesFromLinearState(solveAnalyticIvp(sys.A, sys.c, y0, tau'), sys);
            Vs = max(abs(Van(~locked,:)), [], 'all'); Es = max(abs([e22an; e23an]), [], 'all');
            msg = sprintf('%d cycles, rho = %.3f, p = %d', nreps, rho, sys.p);
            tc.verifyEqual(V(:,~locked)', Van(~locked,:), 'AbsTol', 1e-6*Vs, ['slip rate: ' msg]);
            tc.verifyEqual([e22'; e23'], [e22an; e23an], 'AbsTol', 1e-6*Es, ['strain rate: ' msg]);
        end
    end
end
