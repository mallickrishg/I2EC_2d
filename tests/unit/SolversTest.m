classdef (SharedTestFixtures = {I2ECPathFixture}) SolversTest < matlab.unittest.TestCase
    %SOLVERSTEST  Linear (Newtonian) semi-analytical solvers and their helpers.
    %
    % egh_modal, negateGrowingModes, analyticSolveIVP,
    % computePeriodicInitialConditionBalanced, plus the older
    % computePeriodicInitialCondition / reconditionODEoperator pair,
    % assembleStressKernel, extractComponentsSolutionVector,
    % computeInterseismicSteadystate and computeKLmodesMesh.
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
        %% ---------------- egh_modal ----------------
        function modalIntegralsMatchClosedForm(tc)
            lam = [-1e-9; -3e-10 + 2e-10i; 2e-11; -1e-12];
            t = 3.15e9;
            [e, g, h] = egh_modal(lam, t);
            tc.verifyEqual(e, exp(lam*t), 'RelTol', 1e-12);
            tc.verifyEqual(g, (exp(lam*t) - 1)./lam, 'RelTol', 1e-9);
            tc.verifyEqual(h, ((exp(lam*t) - 1)./lam - t)./lam, 'RelTol', 1e-6);
        end

        function modalIntegralsAreContinuousAtSeriesSwitch(tc)
            t = 1;
            lam = 1e-2 * [1 - 1e-12; 1 + 1e-12; -(1 - 1e-12); -(1 + 1e-12)];
            [~, g, h] = egh_modal(lam, t);
            tc.verifyEqual(g(1), g(2), 'RelTol', 1e-12);
            tc.verifyEqual(h(1), h(2), 'RelTol', 1e-12);
            tc.verifyEqual(g(3), g(4), 'RelTol', 1e-12);
            tc.verifyEqual(h(3), h(4), 'RelTol', 1e-12);
        end

        function modalIntegralsHaveIntegratorLimit(tc)
            [e, g, h] = egh_modal(0, 7);
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

        %% ---------------- analyticSolveIVP ----------------
        function ivpMatchesMatrixExponential(tc)
            y0 = randn(6,1);
            tvec = [0.1 1 10 100 1000] * 3.15e7;
            [Y, U] = analyticSolveIVP(tc.A, tc.c, y0, tvec);
            [Yr, Ur] = SolversTest.expmIVP(tc.A, tc.c, y0, tvec);
            tc.verifyEqual(Y, Yr, 'RelTol', 1e-9, 'AbsTol', 1e-14);
            tc.verifyEqual(U, Ur, 'RelTol', 1e-9, 'AbsTol', 1e-6);
        end

        function ivpAtTimeZeroReturnsInitialState(tc)
            y0 = randn(6,1);
            [Y, U] = analyticSolveIVP(tc.A, tc.c, y0, 0);
            tc.verifyEqual(Y, y0, 'AbsTol', 1e-12);
            tc.verifyEqual(U, zeros(6,1), 'AbsTol', 1e-12);
        end

        %% ---------------- periodic initial condition ----------------
        function periodicStateClosesTheCycle(tc)
            % y(T) + b = y0, checked with the IVP solver
            y0 = computePeriodicInitialConditionBalanced(tc.A, tc.c, tc.b, tc.T);
            yT = analyticSolveIVP(tc.A, tc.c, y0, tc.T);
            tc.verifyEqual(yT + tc.b, y0, 'RelTol', 1e-9, 'AbsTol', 1e-12);
        end

        function balancedAndLegacyPeriodicSolversAgree(tc)
            % for a stable, well-conditioned operator the old expm/SVD
            % solver and the modal solver must give the same state
            y1 = computePeriodicInitialConditionBalanced(tc.A, tc.c, tc.b, tc.T);
            [~, y0] = evalc('computePeriodicInitialCondition(tc.A, tc.c, tc.b, tc.T)');
            tc.verifyEqual(y0, y1, 'RelTol', 1e-8, 'AbsTol', 1e-12);
        end

        function reconditionedOperatorIsStable(tc)
            Ar = reconditionODEoperator(tc.Aunst);
            tc.verifyTrue(isreal(Ar));
            tc.verifyLessThan(max(real(eig(Ar))), 0);
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
            [V2, a, b] = extractComponentsSolutionVector(sol, locked, N, M);
            tc.verifyEqual(V2, V);
            tc.verifyEqual(a, e22);
            tc.verifyEqual(b, e23);
        end

        %% ---------------- computeKLmodesMesh (unused) ----------------
        function klModesAreOrthonormalAndSorted(tc)
            mesh.xc = [linspace(0, 1, 30)' sin(linspace(0, 3, 30))'];
            [V, Vx, Vy, lam] = computeKLmodesMesh(mesh, 0.3);
            tc.verifyEqual(V'*V, eye(size(V,2)), 'AbsTol', 1e-10);
            tc.verifyTrue(issorted(lam, 'descend'));
            tc.verifySize(Vx, size(V)); tc.verifySize(Vy, size(V));
        end

        function klModeDerivativesMatchNystromExtension(tc)
            % Vx(:,k) is the x-derivative of v_k(x) = (1/lam_k) sum_j k(x,x_j) v_k(x_j)
            mesh.xc = [linspace(0, 1, 25)' 0.5*cos(linspace(0, 2, 25))'];
            L = 0.4; m = 4;
            [V, Vx, ~, lam] = computeKLmodesMesh(mesh, L, m);
            h = 1e-6;
            for k = 1:m
                ext = @(p) exp(-((p(1) - mesh.xc(:,1)).^2 + (p(2) - mesh.xc(:,2)).^2) / L^2)' * V(:,k) / lam(k);
                i = 12;
                d = (ext(mesh.xc(i,:) + [h 0]) - ext(mesh.xc(i,:) - [h 0])) / (2*h);
                tc.verifyEqual(Vx(i,k), d, 'RelTol', 1e-5, 'AbsTol', 1e-8);
            end
        end
    end

    methods (Test, TestTags = {'KnownIssue'})
        function reconditionReplacesGrowingModesByMinusEps(tc)
            % help text: eigenvalues with Re >= 0 are replaced by -eps;
            % the code uses -0.1 s^-1 (a 10-second relaxation time in a
            % model with year-to-millennium time scales)
            Ar = reconditionODEoperator(tc.Aunst);
            lam = eig(Ar);
            % every original rate is slower than 1/(3 yr); a replaced mode
            % at -eps is slower still.  None should relax in seconds.
            tc.verifyEqual(sum(real(lam) < -1e-6), 0, ...
                sprintf('min Re(lambda) = %.3g s^-1; documented value is -eps', min(real(lam))));
        end

        function legacyPeriodicSolverIsQuiet(tc)
            out = evalc('computePeriodicInitialCondition(tc.A, tc.c, tc.b, tc.T);');
            tc.verifyEmpty(strtrim(out), 'prints the condition number on every call');
        end

        function interseismicSteadyStateRuns(tc)
            % fails: calls getStrainratesLongterm without the exponent n
            ex = loadExampleModel(true);
            ex.rcv.pinnedPosition = abs(ex.rcv.xc(:,2)) < 30e3;
            f = @() computeInterseismicSteadystate(ex.rcv, ex.shz, ex.evl, 1e-9, 10);
            tc.verifyWarningFree(f);
        end
    end

    methods (Static)
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
