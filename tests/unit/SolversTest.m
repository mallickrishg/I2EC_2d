classdef (SharedTestFixtures = {I2ECPathFixture}) SolversTest < matlab.unittest.TestCase
    %SOLVERSTEST  Linear (Newtonian) semi-analytical solvers and their helpers.
    %
    % computeModalIntegrals, negateGrowingModes, solveAnalyticIvp,
    % computePeriodicInitialCondition, assembleStressKernel and
    % extractRatesFromLinearState.
    %
    % The matrix exponential of the augmented system [A c; 0 0] is the
    % reference solution throughout.  tests/test_negate_exact_solver.m
    % runs the same solvers on the full example operator.

    properties
        A       % 6 x 6 stable, non-normal test operator (1/s)
        Aunst   % same eigenvectors, two eigenvalues flipped to Re > 0
        c
        b
        T = 100*3.15e7
    end

    methods (TestClassSetup)
        function buildOperators(tc)
            rng(7);
            V = randn(6) + 0.3*eye(6);
            yr = 3.15e7;
            lam = -[1/10; 1/50; 1/200; 1/1000; 1/3; 1/30] / yr;
            tc.A = real(V * diag(lam) / V);
            lamU = lam; lamU([2 5]) = -lamU([2 5]);
            tc.Aunst = real(V * diag(lamU) / V);
            tc.c = randn(6,1) * 1e-12;
            tc.b = randn(6,1) * 1e-3;
        end
    end

    methods (Test)
        %% ---------------- computeModalIntegrals ----------------
        function modalIntegralsMatchClosedForm(tc)
            lam = [-1e-9; -3e-10 + 2e-10i; 2e-11; -1e-12];
            t = 3.15e9;
            [e, g, h] = computeModalIntegrals(lam, t);
            tc.verifyEqual(e, exp(lam*t), 'RelTol', 1e-12);
            tc.verifyEqual(g, (exp(lam*t) - 1)./lam, 'RelTol', 1e-9);
            tc.verifyEqual(h, ((exp(lam*t) - 1)./lam - t)./lam, 'RelTol', 1e-6);
        end

        function modalIntegralsAreContinuousAtSeriesSwitch(tc)
            t = 1;
            lam = 1e-2 * [1 - 1e-12; 1 + 1e-12; -(1 - 1e-12); -(1 + 1e-12)];
            [~, g, h] = computeModalIntegrals(lam, t);
            tc.verifyEqual(g(1), g(2), 'RelTol', 1e-12);
            tc.verifyEqual(h(1), h(2), 'RelTol', 1e-12);
            tc.verifyEqual(g(3), g(4), 'RelTol', 1e-12);
            tc.verifyEqual(h(3), h(4), 'RelTol', 1e-12);
        end

        function modalIntegralsHaveIntegratorLimit(tc)
            [e, g, h] = computeModalIntegrals(0, 7);
            tc.verifyEqual([e g h], [1 7 49/2], 'AbsTol', 1e-14);
        end

        %% ---------------- negateGrowingModes ----------------
        function negationKeepsStableOperator(tc)
            [~, lam, An] = negateGrowingModes(tc.A);
            tc.verifyLessThan(max(real(lam)), 0);
            tc.verifyEqual(An, tc.A, 'AbsTol', 1e-8*norm(tc.A));
        end

        function negationFlipsGrowingModesOnly(tc)
            [~, lam, An] = negateGrowingModes(tc.Aunst);
            tc.verifyLessThan(max(real(lam)), 0);
            tc.verifyTrue(isreal(An));
            tc.verifyEqual(sort(real(eig(An))), sort(real(eig(tc.A))), 'RelTol', 1e-8);
        end

        function solversRejectGrowingModes(tc)
            % the solvers no longer negate: an operator with growing modes
            % (centroid kernels) must be negated explicitly first
            y0 = ones(size(tc.c));
            tc.verifyError(@() solveAnalyticIvp(tc.Aunst, tc.c, y0, tc.T), 'I2EC:growingModes');
            tc.verifyError(@() computePeriodicInitialCondition(tc.Aunst, tc.c, tc.b, tc.T), 'I2EC:growingModes');
            [~, ~, An] = negateGrowingModes(tc.Aunst);
            tc.verifyWarningFree(@() solveAnalyticIvp(An, tc.c, y0, tc.T));
            tc.verifyWarningFree(@() computePeriodicInitialCondition(An, tc.c, tc.b, tc.T));
        end

        function periodicSolverRejectsNeutralModes(tc)
            % a mode that neither decays nor grows (an integrator, or an
            % undamped oscillation) leaves the periodic state undetermined
            A0 = blkdiag(tc.A, 0);
            tc.verifyError(@() computePeriodicInitialCondition(A0, [tc.c; 0], [tc.b; 0], tc.T), ...
                'I2EC:neutralModes');
            w = 2*pi/(7*tc.T);
            Aosc = blkdiag(tc.A, [0 w; -w 0]);
            tc.verifyError(@() computePeriodicInitialCondition(Aosc, [tc.c; 0; 0], [tc.b; 0; 0], tc.T), ...
                'I2EC:neutralModes');
            % the IVP is still defined
            tc.verifyWarningFree(@() solveAnalyticIvp(A0, [tc.c; 1e-12], [ones(size(tc.c)); 0], tc.T));
        end

        function periodicSolverAcceptsWidelySeparatedTimescales(tc)
            % a strong mantle with a weak fault: modes 1e-11 of max|lambda|
            % that still relax over many cycles (|lambda|*T ~ 1e-9) are
            % solved, not rejected as neutral, and close the cycle
            As = blkdiag(tc.A, 1e-8*tc.A);
            cs = [tc.c; 1e-8*tc.c]; bs = [tc.b; tc.b];
            y0 = computePeriodicInitialCondition(As, cs, bs, tc.T);
            yT = solveAnalyticIvp(As, cs, y0, tc.T);
            tc.verifyLessThan(norm(yT + bs - y0) / norm(y0), 1e-6);
        end

        function exampleOperatorIsAccurateForWeakAndStrongMantles(tc)
            % The example 2 operator with a moderate (t_M/T = 0.1) and a
            % strong (1e8) mantle; no warning from the periodic solve. Fault and mantle rates differ by ~1e7 in
            % size, so each eigendecomposition alone loses digits in one of
            % the two regimes (unbalanced: mantle 1e-9 at 0.1; balanced:
            % fault 2e-9 for strong mantles). The periodic state is solved
            % directly and the evolution keeps the more accurate modes.
            ex = loadExampleModel(true);
            yr = 3.15e7; T = 500*yr;
            for ratio = [0.1 1e8]
                [A, c, b, Nc] = SolversTest.exampleSystem(ex, ratio*T*30e3, T);
                n = size(A, 1); f = 1:Nc; mt = Nc+1:n;
                y0 = tc.verifyWarningFree(@() computePeriodicInitialCondition(A, c, b, T));
                % one cycle from y0 returns to y0 - b
                yT = solveAnalyticIvp(A, c, y0, T);
                tc.verifyLessThan(norm(yT(f) + b(f) - y0(f)) / norm(y0(f)), 1e-10, sprintf('fault, t_M/T = %g', ratio));
                tc.verifyLessThan(norm(yT(mt) + b(mt) - y0(mt)) / norm(y0(mt)), 1e-10, sprintf('mantle, t_M/T = %g', ratio));
                % evolution against the matrix exponential
                t = 1e-2*T;
                X = expm([A c; zeros(1, n+1)] * t);
                yref = X(1:n, 1:n)*y0 + X(1:n, end);
                y = solveAnalyticIvp(A, c, y0, t);
                tc.verifyLessThan(norm(y(f) - yref(f)) / norm(yref(f)), 1e-11, sprintf('evolution, fault, t_M/T = %g', ratio));
                tc.verifyLessThan(norm(y(mt) - yref(mt)) / norm(yref(mt)), 1e-11, sprintf('evolution, mantle, t_M/T = %g', ratio));
            end
        end

        function eigOperatorOptionalBalancingGivesSameModes(tc)
            % unbalanced by default; balancing is an option with the same
            % eigenvalues and eigenvectors of A
            [V0, l0] = eigOperator(tc.A);
            [V1, l1] = eigOperator(tc.A, [], true);
            tc.verifyEqual(sort(l1, 'ComparisonMethod', 'real'), sort(l0, 'ComparisonMethod', 'real'), 'RelTol', 1e-10);
            tc.verifyEqual(real(V0*diag(l0)/V0), tc.A, 'AbsTol', 1e-10*norm(tc.A));
            tc.verifyEqual(real(V1*diag(l1)/V1), tc.A, 'AbsTol', 1e-10*norm(tc.A));
        end

        function eigOperatorMatchesEig(tc)
            [V, lam, ng] = eigOperator(tc.A);
            tc.verifyEqual(sort(lam, 'ComparisonMethod', 'real'), sort(eig(tc.A), 'ComparisonMethod', 'real'), 'RelTol', 1e-8);
            tc.verifyEqual(real(V*diag(lam)/V), tc.A, 'AbsTol', 1e-8*norm(tc.A));
            tc.verifyEqual(ng, 0);
            [~, ~, ng] = eigOperator(tc.Aunst);
            tc.verifyEqual(ng, 2);
        end

        %% ---------------- solveAnalyticIvp ----------------
        function ivpMatchesMatrixExponential(tc)
            y0 = randn(6,1);
            tvec = [0.1 1 10 100 1000] * 3.15e7;
            [Y, U] = solveAnalyticIvp(tc.A, tc.c, y0, tvec);
            [Yr, Ur] = SolversTest.expmIVP(tc.A, tc.c, y0, tvec);
            tc.verifyEqual(Y, Yr, 'RelTol', 1e-9, 'AbsTol', 1e-14);
            tc.verifyEqual(U, Ur, 'RelTol', 1e-9, 'AbsTol', 1e-6);
        end

        function ivpAtTimeZeroReturnsInitialState(tc)
            y0 = randn(6,1);
            [Y, U] = solveAnalyticIvp(tc.A, tc.c, y0, 0);
            tc.verifyEqual(Y, y0, 'AbsTol', 1e-12);
            tc.verifyEqual(U, zeros(6,1), 'AbsTol', 1e-12);
        end

        %% ---------------- periodic initial condition ----------------
        function periodicStateClosesTheCycle(tc)
            % y(T) + b = y0, checked with the IVP solver
            y0 = computePeriodicInitialCondition(tc.A, tc.c, tc.b, tc.T);
            yT = solveAnalyticIvp(tc.A, tc.c, y0, tc.T);
            tc.verifyEqual(yT + tc.b, y0, 'RelTol', 1e-9, 'AbsTol', 1e-12);
        end

        %% ---------------- bookkeeping ----------------
        function assembledKernelHasExpectedBlocks(tc)
            ex = loadExampleModel(true);
            locked = abs(ex.rcv.xc(:,2)) < 30e3;
            K = assembleStressKernel(ex.evl, locked);
            nf = sum(~locked); ns = ex.shz.N;
            tc.verifySize(K, [nf + 2*ns, nf + 2*ns]);
            tc.verifyEqual(K(1:nf, 1:nf), ex.evl.KK(~locked, ~locked));
            tc.verifyEqual(K(nf+1:nf+ns, nf+ns+1:end), ex.evl.LL(:,:,1,2));
            tc.verifyEqual(K(nf+ns+1:end, 1:nf), ex.evl.KL(:, ~locked, 2));
        end

        function solutionVectorRoundTrip(tc)
            N = 8; M = 5; Nt = 3;
            locked = false(N,1); locked([1 2 7]) = true;
            V = rand(N, Nt); V(locked,:) = 0;
            e22 = rand(M, Nt); e23 = rand(M, Nt);
            sol = [V(~locked,:); e22; e23];
            [V2, a, b] = extractRatesFromLinearState(sol, locked, N, M);
            tc.verifyEqual(V2, V);
            tc.verifyEqual(a, e22);
            tc.verifyEqual(b, e23);
        end
    end

    methods (Static)
        function [A, c, b, Nc] = exampleSystem(ex, eta, T)
            % example 2's linear system with a uniform mantle viscosity eta (MPa s)
            rcv = ex.rcv; shz = ex.shz; evl = ex.evl; Vpl = 1e-9;
            rcv.faultLaw = 'linear';
            rcv.Asigma = 1e-6*(5e17/sum(rcv.W)) .* ones(rcv.N, 1);
            locked = abs(rcv.xc(:,2)) < 30e3;
            rcv.pinnedPosition = locked; rcv.Vpl = Vpl .* ones(rcv.N, 1);
            shz.n = ones(shz.N, 1); shz.alpha = (1/eta) .* ones(shz.N, 1);
            [e22, e23] = computeLongTermStrainRates(shz, rcv.dip(1)*pi/180, [0, 20e3], [-140e3, 35e3], 1);
            shz.e22pl = e22 .* Vpl; shz.e23pl = -e23 .* Vpl;
            slip = T*Vpl .* locked;
            sc = struct('Nevents', 1, 'Timing', 0, 'dtau', evl.KK*slip, ...
                'dsigma22', evl.KL(:,:,1)*slip, 'dsigma23', evl.KL(:,:,2)*slip);
            sys = assembleLinearSystem(evl, rcv, shz, locked, sc);
            A = sys.A; c = sys.c; b = sys.b; Nc = sys.Nc;
        end

        function [Y, U] = expmIVP(A, c, y0, tvec)
            % reference: y and its time integral from expm of the augmented system
            n = numel(y0);
            M = zeros(2*n + 1);
            M(1:n, 1:n) = A;
            M(1:n, 2*n+1) = c;
            M(n+1:2*n, 1:n) = eye(n);
            Y = zeros(n, numel(tvec)); U = Y;
            for j = 1:numel(tvec)
                z = expm(M*tvec(j)) * [y0; zeros(n,1); 1];
                Y(:,j) = z(1:n);
                U(:,j) = z(n+1:2*n);
            end
        end
    end
end
