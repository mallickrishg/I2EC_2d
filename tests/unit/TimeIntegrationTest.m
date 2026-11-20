classdef (SharedTestFixtures = {I2ECPathFixture}) TimeIntegrationTest < matlab.unittest.TestCase
    %TIMEINTEGRATIONTEST  Numerical earthquake-cycle path (MAIN_i2ec_*):
    % ode_viscoelastic, viscoeq_sequence_simulation,
    % runImposedEarthquakeCycles, extract_v_edot_from_statevector.
    %
    % Key check: with the whole fault locked and Newtonian mantle the
    % numerical path is linear, so ode45 must reproduce the semi-analytical
    % solver (analyticSolveIVP) through a coseismic step.  This ties the
    % two solution routes of the code together.

    properties
        M       % small model
        Evl     % kernels with kernelmodify = 1 (no growing modes)
        yr = 3.15e7
        Vpl = 1e-9
    end

    methods (TestClassSetup)
        function build(tc)
            tc.M = makeSmallModel(tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder);
            [~, tc.Evl] = evalc('computeAllStressKernelsBem(tc.M.rcv, tc.M.shz, tc.M.bnd, ''kernelmodify'', 1)');
        end
    end

    methods (TestMethodSetup)
        function setRheology(tc)
            rcv = tc.M.rcv; shz = tc.M.shz;
            rcv.Vpl = tc.Vpl * ones(rcv.N, 1);
            rcv.Asigma = 0.5 * ones(rcv.N, 1);
            rcv.pinnedPosition = false(rcv.N, 1);
            shz.n = ones(shz.N, 1);
            shz.alpha = 1/(1e19*1e-6) * ones(shz.N, 1);
            shz.e22pl = 1e-15 * linspace(-1, 1, shz.N)';
            shz.e23pl = 2e-15 * ones(shz.N, 1);
        end
    end

    methods (Test)
        function stateVectorRoundTrip(tc)
            rcv = tc.M.rcv; shz = tc.M.shz;
            V = tc.Vpl * (1:rcv.N); s22 = linspace(1, 2, shz.N); s23 = linspace(-1, 0, shz.N);
            Y = zeros(1, rcv.N + 2*shz.N);
            Y(1:rcv.N) = log(V);
            Y(rcv.N+1:2:end) = s22;
            Y(rcv.N+2:2:end) = s23;
            [V2, e22, e23] = extract_v_edot_from_statevector(Y, rcv, shz);
            tc.verifyEqual(V2, V, 'RelTol', 1e-12);
            tc.verifyEqual(e22, shz.alpha' .* s22, 'RelTol', 1e-12);
            tc.verifyEqual(e23, shz.alpha' .* s23, 'RelTol', 1e-12);
        end

        function steadyStateHasZeroRate(tc)
            % plate-rate slip and long-term strain rate: nothing changes
            rcv = tc.M.rcv; shz = tc.M.shz;
            Y = [log(rcv.Vpl); reshape([shz.e22pl shz.e23pl]' ./ shz.alpha', [], 1)];
            Yp = ode_viscoelastic(0, Y, rcv, shz, tc.Evl);
            % stress-rate terms are ~1e-11 MPa/s each; cancellation leaves round-off
            tc.verifyEqual(Yp, zeros(size(Y)), 'AbsTol', 1e-20);
        end

        function steadyStateSurvivesACycleWithoutEarthquake(tc)
            rcv = tc.M.rcv; shz = tc.M.shz;
            sc = struct('Nevents', 1, 'Timing', 10*tc.yr, ...
                'dtau', zeros(rcv.N,1), 'dsigma22', zeros(shz.N,1), 'dsigma23', zeros(shz.N,1));
            [~, t, V, e22, e23] = evalc('runImposedEarthquakeCycles(rcv, shz, tc.Evl, sc, 2, 50*tc.yr)');
            tc.verifyEqual(t(end), 50*tc.yr, 'RelTol', 1e-12);
            tc.verifyEqual(V, repmat(rcv.Vpl', size(V,1), 1), 'RelTol', 1e-8);
            tc.verifyEqual(e22, repmat(shz.e22pl', size(V,1), 1), 'AbsTol', 1e-22);
            tc.verifyEqual(e23, repmat(shz.e23pl', size(V,1), 1), 'AbsTol', 1e-22);
        end

        function eventSequenceCoversTheCycle(tc)
            rcv = tc.M.rcv; shz = tc.M.shz;
            Trecur = 100*tc.yr;
            sc = struct('Nevents', 2, 'Timing', [5 40]*tc.yr, ...
                'dtau', zeros(rcv.N,2), 'dsigma22', zeros(shz.N,2), 'dsigma23', zeros(shz.N,2));
            Y0 = [log(rcv.Vpl); reshape([shz.e22pl shz.e23pl]' ./ shz.alpha', [], 1)];
            [t, Y] = viscoeq_sequence_simulation(@(t,y) ode_viscoelastic(t,y,rcv,shz,tc.Evl), Y0, rcv, shz, sc, Trecur);
            tc.verifyEqual(t(1), 0);
            tc.verifyEqual(t(end), Trecur, 'RelTol', 1e-12);
            tc.verifyTrue(all(diff(t) >= 0), 'time is non-decreasing');
            tc.verifyEqual(size(Y,1), numel(t));
        end

        function cyclesChainFromTheEndState(tc)
            % nreps = 3 returns the cycle that starts from the end state of
            % nreps = 2: same ode45 calls in the same order, so bitwise equal
            [rcv, shz, sc, Y0, f, Trecur] = tc.multiEventSetup();
            [~, ~, Y2] = evalc('viscoeq_sequence_simulation(f, Y0, rcv, shz, sc, Trecur, 2)');
            [~, t3, Y3] = evalc('viscoeq_sequence_simulation(f, Y0, rcv, shz, sc, Trecur, 3)');
            [t1, Y1] = viscoeq_sequence_simulation(f, Y2(end,:)', rcv, shz, sc, Trecur);
            tc.verifyEqual(t3, t1);
            tc.verifyEqual(Y3, Y1);
            tc.verifyNotEqual(Y2(end,:), Y0', 'the events do change the state');
        end

        function everyEventIsAppliedOnce(tc)
            % K events per cycle: each event time appears twice in t (just
            % before and after the step), and the jump in log(v) and in the
            % mantle stress is exactly that event's stress change
            [rcv, shz, sc, Y0, f, Trecur] = tc.multiEventSetup();
            [t, Y] = viscoeq_sequence_simulation(f, Y0, rcv, shz, sc, Trecur);
            iV = 1:rcv.N; i22 = rcv.N+1:2:size(Y,2); i23 = rcv.N+2:2:size(Y,2);
            for k = 1:sc.Nevents
                i = find(t == sc.Timing(k));
                tc.verifyNumElements(i, 2, sprintf('event %d', k));
                dY = Y(i(2),:) - Y(i(1),:);
                tc.verifyEqual(dY(iV)', sc.dtau(:,k)./rcv.Asigma, 'AbsTol', 1e-12);
                tc.verifyEqual(dY(i22)', sc.dsigma22(:,k), 'AbsTol', 1e-12);
                tc.verifyEqual(dY(i23)', sc.dsigma23(:,k), 'AbsTol', 1e-12);
            end
            tc.verifyEqual(t(end), Trecur, 'RelTol', 1e-12);
            tc.verifyTrue(all(diff(t) >= 0));
        end

        function badEventTimesAreRejected(tc)
            [rcv, shz, sc, Y0, f, Trecur] = tc.multiEventSetup();
            run = @(T) viscoeq_sequence_simulation(f, Y0, rcv, shz, ...
                setfield(sc, 'Timing', T), Trecur); %#ok<SFLD>
            tc.verifyError(@() run(sc.Timing([2 1 3 4 5])), 'viscoeq_sequence_simulation:timing');
            tc.verifyError(@() run([0 sc.Timing(2:end)]), 'viscoeq_sequence_simulation:timing');
            tc.verifyError(@() run([sc.Timing(1:end-1) Trecur]), 'viscoeq_sequence_simulation:timing');
            tc.verifyError(@() run(sc.Timing(1:4)), 'viscoeq_sequence_simulation:events');
        end

        function lockedNewtonianCycleMatchesAnalyticSolver(tc)
            % fully locked fault + n = 1 mantle: the ODE is linear in the
            % mantle strain rate e = alpha * s:
            %   de/dt = alpha .* (LL (e - e_pl) - KL V_pl) = A e + c
            rcv = tc.M.rcv; shz = tc.M.shz; evl = tc.Evl;
            rcv.pinnedPosition = true(rcv.N, 1);
            Ns = shz.N; alpha = [shz.alpha; shz.alpha];
            LL = [evl.LL(:,:,1,1) evl.LL(:,:,1,2); evl.LL(:,:,2,1) evl.LL(:,:,2,2)];
            KL = [evl.KL(:,:,1); evl.KL(:,:,2)];
            epl = [shz.e22pl; shz.e23pl];
            A = alpha .* LL;
            c = -A*epl - alpha .* (KL*rcv.Vpl);

            % one event at t1 with coseismic slip of 2 m on the whole fault
            t1 = 10*tc.yr; Trecur = 60*tc.yr;
            slip = 2 * ones(rcv.N, 1);
            sc = struct('Nevents', 1, 'Timing', t1, 'dtau', zeros(rcv.N,1), ...
                'dsigma22', evl.KL(:,:,1)*slip, 'dsigma23', evl.KL(:,:,2)*slip);
            e0 = epl;
            Y0 = [log(rcv.Vpl); reshape([e0(1:Ns) e0(Ns+1:end)]' ./ shz.alpha', [], 1)];
            [t, Y] = viscoeq_sequence_simulation(@(t,y) ode_viscoelastic(t,y,rcv,shz,evl), Y0, rcv, shz, sc, Trecur);
            [~, e22, e23] = extract_v_edot_from_statevector(Y(end,:), rcv, shz);

            % analytic: evolve to t1, add coseismic jump, evolve to Trecur
            e1 = analyticSolveIVP(A, c, e0, t1);
            e1 = e1 + alpha .* (KL*slip);
            eT = analyticSolveIVP(A, c, e1, Trecur - t1);
            tc.verifyEqual(t(end), Trecur, 'RelTol', 1e-12);
            tc.verifyEqual([e22'; e23'], eT, 'RelTol', 1e-5, 'AbsTol', 1e-6*max(abs(eT)));
        end

        function spunUpLockedNewtonianCycleMatchesAnalyticPeriodicCycle(tc)
            % Numerical path with spin-up (runImposedEarthquakeCycles,
            % ode45) against the semi-analytical periodic cycle
            % (computePeriodicInitialConditionBalanced + analyticSolveIVP).
            % Fully locked fault + n = 1 mantle, so both paths solve the
            % same linear ODE de/dt = A e + c (e = alpha s). Both use the
            % same operator: the unmodified LL with its growing modes
            % negated (negateGrowingModes), so every mode decays and the
            % spin-up converges to the periodic cycle at the rate of the
            % slowest mode, rho = max|exp(lambda T)| per cycle.
            rcv = tc.M.rcv; shz = tc.M.shz;
            rcv.pinnedPosition = true(rcv.N, 1);
            shz.alpha = 1/(1e18*1e-6) * ones(shz.N, 1);   % slowest mode ~36 yr
            Ns = shz.N; alpha = [shz.alpha; shz.alpha];
            [~, evl] = evalc('computeAllStressKernelsBem(tc.M.rcv, tc.M.shz, tc.M.bnd, ''kernelmodify'', 0)');
            LL = [evl.LL(:,:,1,1) evl.LL(:,:,1,2); evl.LL(:,:,2,1) evl.LL(:,:,2,2)];
            [~, lambda, A] = negateGrowingModes(alpha .* LL);
            LLneg = A ./ alpha;                            % same operator for ode45
            evl.LL(:,:,1,1) = LLneg(1:Ns, 1:Ns);     evl.LL(:,:,1,2) = LLneg(1:Ns, Ns+1:end);
            evl.LL(:,:,2,1) = LLneg(Ns+1:end, 1:Ns); evl.LL(:,:,2,2) = LLneg(Ns+1:end, Ns+1:end);
            KL = [evl.KL(:,:,1); evl.KL(:,:,2)];
            epl = [shz.e22pl; shz.e23pl];
            c = -A*epl - alpha .* (KL*rcv.Vpl);

            % one event per cycle at t1: 2 m of slip on the whole fault
            Trecur = 100*tc.yr; t1 = 30*tc.yr;
            slip = 2 * ones(rcv.N, 1);
            sc = struct('Nevents', 1, 'Timing', t1, 'dtau', zeros(rcv.N,1), ...
                'dsigma22', evl.KL(:,:,1)*slip, 'dsigma23', evl.KL(:,:,2)*slip);
            b = alpha .* (KL*slip);

            % analytic periodic cycle: y0 just after the event, y(T) + b = y0
            y0 = computePeriodicInitialConditionBalanced(A, c, b, Trecur);

            % enough cycles for the start-up transient to fall below 1e-8
            rho = max(abs(exp(lambda*Trecur)));
            nreps = ceil(log(1e-8)/log(rho)) + 1;
            tc.assertLessThan(nreps, 20, sprintf('spin-up too slow: rho = %.3f', rho));
            [~, t, ~, e22, e23] = evalc('runImposedEarthquakeCycles(rcv, shz, evl, sc, nreps, Trecur)');
            Enum = [e22'; e23'];

            % analytic state at each output time: time since the last event
            post = [false; diff(t) == 0];                 % post-event copy of t1
            tau = t - t1 + Trecur*(t < t1 | (t == t1 & ~post));
            Ean = analyticSolveIVP(A, c, y0, tau');
            scale = max(abs(Ean(:)));
            tc.verifyEqual(Enum, Ean, 'AbsTol', 1e-6*scale, ...
                sprintf('spun-up cycle (%d cycles, rho = %.3f) vs analytic periodic cycle', nreps, rho));

            % without spin-up the first cycle is still far from periodic
            % (compared at the end of the cycle)
            [~, ~, ~, f22, f23] = evalc('runImposedEarthquakeCycles(rcv, shz, evl, sc, 1, Trecur)');
            gap1 = max(abs([f22(end,:)'; f23(end,:)'] - Ean(:,end)));
            tc.verifyGreaterThan(gap1, 1e-3*scale, 'one cycle from the long-term state is not yet periodic');
        end
    end

    methods (Test)
        function singleCycleRunReturnsTime(tc)
            % nreps = 1: one cycle from the long-term state, no spin-up
            % (the first call used to discard t, so t was never assigned)
            rcv = tc.M.rcv; shz = tc.M.shz;
            sc = struct('Nevents', 1, 'Timing', 10*tc.yr, ...
                'dtau', zeros(rcv.N,1), 'dsigma22', zeros(shz.N,1), 'dsigma23', zeros(shz.N,1));
            Trecur = 50*tc.yr;
            [~, t, V, e22] = evalc('runImposedEarthquakeCycles(rcv, shz, tc.Evl, sc, 1, Trecur)');
            tc.verifyEqual(t(end), Trecur, 'RelTol', 1e-12);
            tc.verifySize(V, [numel(t) rcv.N]);
            tc.verifySize(e22, [numel(t) shz.N]);
        end
    end

    methods (Access = private)
        function [rcv, shz, sc, Y0, f, Trecur] = multiEventSetup(tc)
            % 5 irregularly spaced events with different slip, n = 3
            % mantle, a third of the fault locked
            rcv = tc.M.rcv; shz = tc.M.shz; evl = tc.Evl;
            shz.n = 3 * ones(shz.N, 1);
            rcv.pinnedPosition = (1:rcv.N)' <= rcv.N/3;
            Trecur = 80*tc.yr;
            rng(7); slip = rand(rcv.N, 5);
            sc = struct('Nevents', 5, 'Timing', [3 11 12 40 66]*tc.yr, ...
                'dtau', 0.3*rand(rcv.N, 5), ...
                'dsigma22', evl.KL(:,:,1)*slip, 'dsigma23', evl.KL(:,:,2)*slip);
            Y0 = [log(rcv.Vpl); reshape([shz.e22pl shz.e23pl]' ./ shz.alpha', [], 1)];
            f = @(t, y) ode_viscoelastic(t, y, rcv, shz, evl);
        end
    end
end
