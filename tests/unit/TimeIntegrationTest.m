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
    end

    methods (Test, TestTags = {'KnownIssue'})
        function singleCycleRunReturnsTime(tc)
            % with nreps = 1 the first call discards t ([~,Y] = ...) and the
            % loop never runs, so output t is never assigned
            rcv = tc.M.rcv; shz = tc.M.shz;
            sc = struct('Nevents', 1, 'Timing', 10*tc.yr, ...
                'dtau', zeros(rcv.N,1), 'dsigma22', zeros(shz.N,1), 'dsigma23', zeros(shz.N,1));
            msg = '';
            try
                evalc('[t, V] = runImposedEarthquakeCycles(rcv, shz, tc.Evl, sc, 1, 50*tc.yr);');
            catch err
                msg = err.message;
            end
            tc.verifyEmpty(msg, msg);
        end
    end
end
