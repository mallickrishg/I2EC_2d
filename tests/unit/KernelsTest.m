classdef (SharedTestFixtures = {I2ECPathFixture}) KernelsTest < matlab.unittest.TestCase
    %KERNELSTEST  BEM stress and displacement kernels (computeAll*KernelsBem).
    %
    % Uses the 10-patch / 24-triangle synthetic model from makeSmallModel
    % (about a second to build) and the cached example kernels in kernels/.
    %
    % The central check is Betti's reciprocal theorem: with element sizes
    % as weights, the fault-fault, fault-mantle and mantle-mantle blocks
    % must be (approximately, for well-separated elements) symmetric.  A
    % sign or component-ordering error in any block breaks it.

    properties
        Sm     % small model struct (makeSmallModel)
        Evl    % small-model stress kernels, kernelmodify = 0
        Devl   % small-model displacement kernels
        Ex     % example model with cached evl_orig / devl
    end

    methods (TestClassSetup)
        function buildKernels(tc)
            folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            tc.Sm = makeSmallModel(folder);
            [~, tc.Evl] = evalc('computeAllStressKernelsBem(tc.Sm.rcv, tc.Sm.shz, tc.Sm.bnd, ''kernelmodify'', 0)');
            [~, tc.Devl] = evalc('computeAllDisplacementKernelsBem(tc.Sm.obs, tc.Sm.rcv, tc.Sm.shz, tc.Sm.bnd, 1)');
            tc.Ex = loadExampleModel(true);
        end
    end

    methods (Test)
        function stressKernelShapes(tc)
            Nf = tc.Sm.rcv.N; Ns = tc.Sm.shz.N;
            tc.verifySize(tc.Evl.KK, [Nf Nf]);
            tc.verifySize(tc.Evl.KL, [Ns Nf 2]);
            tc.verifySize(tc.Evl.LK, [Nf Ns 2]);
            tc.verifySize(tc.Evl.LL, [Ns Ns 2 2]);
            tc.verifySize(tc.Devl.KO, [size(tc.Sm.obs,1) Nf 2]);
            tc.verifySize(tc.Devl.LO, [size(tc.Sm.obs,1) Ns 2 2]);
        end

        function selfStressOpposesDeformation(tc)
            for e = {tc.Evl, tc.Ex.evl}
                tc.verifyTrue(all(diag(e{1}.KK) < 0), 'KK diagonal');
                tc.verifyTrue(all(diag(e{1}.LL(:,:,1,1)) < 0), 'LL e22 diagonal');
                tc.verifyTrue(all(diag(e{1}.LL(:,:,2,2)) < 0), 'LL e23 diagonal');
            end
        end

        function faultFaultKernelIsReciprocal(tc)
            for c = {{tc.Evl, tc.Sm}, {tc.Ex.evl, tc.Ex}}
                [e, m] = c{1}{:};
                W = m.rcv.W; xc = m.rcv.xc;
                far = pairDistance(xc, xc) > 3*max(W);
                M = e.KK .* W;
                r = reciprocityRatio(M, M.', far);
                tc.verifyEqual(median(r), 1, 'AbsTol', 5e-3);
            end
        end

        function faultMantleKernelsAreReciprocal(tc)
            % W_i LK(i,j,c) = 2 a_j KL(j,i,c): the factor 2 comes from the
            % deviatoric work (s22 - s33) e22 + 2 s23 e23
            for c = {{tc.Evl, tc.Sm}, {tc.Ex.evl, tc.Ex}}
                [e, m] = c{1}{:};
                W = m.rcv.W; a = m.shz.area;
                far = pairDistance(m.rcv.xc, m.shz.xc) > 3*(max(W) + sqrt(max(a)));
                for comp = 1:2
                    lhs = e.LK(:,:,comp) .* W;
                    rhs = 2 * (e.KL(:,:,comp) .* a).';
                    r = reciprocityRatio(lhs, rhs, far);
                    tc.verifyEqual(median(r), 1, 'AbsTol', 1e-2, sprintf('component %d', comp));
                    % centroid collocation: a few % spread even for distant pairs
                    tc.verifyGreaterThan(prctile(r, 5), 0.9, sprintf('component %d', comp));
                end
            end
        end

        function mantleMantleKernelIsReciprocal(tc)
            e = tc.Ex.evl; a = tc.Ex.shz.area;
            far = pairDistance(tc.Ex.shz.xc, tc.Ex.shz.xc) > 5*sqrt(max(a));
            L11 = e.LL(:,:,1,1) .* a; L22 = e.LL(:,:,2,2) .* a;
            L12 = e.LL(:,:,1,2) .* a; L21 = e.LL(:,:,2,1) .* a;
            tc.verifyEqual(median(reciprocityRatio(L11, L11.', far)), 1, 'AbsTol', 2e-2);
            tc.verifyEqual(median(reciprocityRatio(L22, L22.', far)), 1, 'AbsTol', 2e-2);
            tc.verifyEqual(median(reciprocityRatio(L12, L21.', far)), 1, 'AbsTol', 5e-2);
        end

        function siteVelocitiesMatchDisplacementKernels(tc)
            % computeSiteVelocitiesBem (unused) vs devl kernels: same BEM problem
            m = tc.Sm;
            rng(1);
            V = rand(2, m.rcv.N); e22 = randn(2, m.shz.N); e23 = randn(2, m.shz.N);
            gps = computeSiteVelocitiesBem(m.obs, m.rcv, m.shz, m.bnd, V, e22, e23);
            vx = (tc.Devl.KO(:,:,1)*V' + tc.Devl.LO(:,:,1,1)*e22' + tc.Devl.LO(:,:,1,2)*e23')';
            vz = (tc.Devl.KO(:,:,2)*V' + tc.Devl.LO(:,:,2,1)*e22' + tc.Devl.LO(:,:,2,2)*e23')';
            tc.verifyEqual(gps.vx, vx, 'AbsTol', 1e-9*max(abs(vx(:))));
            tc.verifyEqual(gps.vz, vz, 'AbsTol', 1e-9*max(abs(vz(:))));
        end

        function kernelModifyRemovesGrowingModes(tc)
            [~, evm] = evalc('computeAllStressKernelsBem(tc.Sm.rcv, tc.Sm.shz, tc.Sm.bnd, ''kernelmodify'', 1)');
            L = [evm.LL(:,:,1,1) evm.LL(:,:,1,2); evm.LL(:,:,2,1) evm.LL(:,:,2,2)];
            scale = max(abs(eig(L)));
            tc.verifyLessThanOrEqual(max(real(eig(L))) / scale, 1e-10);
            tc.verifyEqual(evm.KK, tc.Evl.KK, 'other blocks unchanged');
            tc.verifyEqual(evm.LK, tc.Evl.LK);
        end

        function unknownOptionIsRejected(tc)
            m = tc.Sm;
            % the option is parsed before any kernel is computed
            tc.verifyError(@() computeAllStressKernelsBem(m.rcv, m.shz, m.bnd, 'kernelmodfy', 0), ?MException);
        end

        function faultTractionKernelsBemNeedsTwoOrThreeOutputs(tc)
            m = tc.Sm;
            tc.verifyError(@() computeFaultTractionKernelsBem(m.rcv, m.rcv, m.bnd), ?MException);
        end

        function cachedKernelsMatchCurrentMesh(tc)
            % kernels/*.mat must be built on the meshes in inputs/
            e = tc.Ex.evl;
            tc.verifySize(e.KK, [tc.Ex.rcv.N tc.Ex.rcv.N]);
            tc.verifySize(e.LL, [tc.Ex.shz.N tc.Ex.shz.N 2 2]);
            tc.verifySize(tc.Ex.devl.KO, [size(tc.Ex.obs,1) tc.Ex.rcv.N 2]);
        end
    end

    methods (Test, TestTags = {'KnownIssue'})
        function unmodifiedMantleKernelHasNoGrowingModes(tc)
            % An elastic relaxation operator should not have eigenvalues
            % with positive real part.  The collocated LL block does (43
            % on the example mesh, even on the 24-triangle mesh); the
            % solvers negate them (negateGrowingModes) and
            % kernelmodify = 1 zeroes them.  Kept as a known issue of the
            % discretisation, not of the solvers.
            e = tc.Ex.evl;
            L = [e.LL(:,:,1,1) e.LL(:,:,1,2); e.LL(:,:,2,1) e.LL(:,:,2,2)];
            lam = eig(L);
            tc.verifyEqual(sum(real(lam) > 0), 0, ...
                sprintf('%d eigenvalues with Re > 0, max Re = %.3g', sum(real(lam) > 0), max(real(lam))));
        end

        function stressKernelsIgnoreBoundaryVelocity(tc)
            % Kernels are per unit source and must be linear in it.
            % computeFaultTractionKernelsBem and computeShzStressKernelsBem
            % put boundary.Vx / Vz into the Dirichlet data of every column,
            % so a boundary with non-zero prescribed velocity (bcx/bcz in
            % the .seg file, e.g. boundarytight2d.seg) adds a constant to
            % each column.  The displacement-kernel routines use 0*Vx.
            % own copy of the model: receiver objects are handles
            m = makeSmallModel(tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder);
            [K0, ~] = computeFaultTractionKernelsBem(m.rcv, m.rcv, m.bnd);
            % non-uniform (a uniform value is a rigid translation: no stress)
            m.bnd.Vx = linspace(0, 1, m.bnd.N)';
            [K1, ~] = computeFaultTractionKernelsBem(m.rcv, m.rcv, m.bnd);
            tc.verifyEqual(K1, K0, 'RelTol', 1e-12);
        end
    end

    methods (Test, TestTags = {'Slow'})
        function cachedStressKernelsAreReproducible(tc)
            % rebuild evl_orig from inputs/ and compare to kernels/evl_orig.mat
            % (minutes; catches a stale cache after a mesh or kernel change)
            m = tc.Ex;
            [~, e] = evalc('computeAllStressKernelsBem(m.rcv, m.shz, m.bnd, ''kernelmodify'', 0)');
            % KK, KL come from closed-form line-dislocation stresses;
            % LK, LL from the finite-difference triangle stress (step 1e-6 of
            % the shortest edge), which reproduces only to ~1e-5 (measured
            % max 1e-5, median 3e-6 relative, 2026-11-19)
            tol = struct('KK', 1e-8, 'KL', 1e-8, 'LK', 1e-4, 'LL', 1e-4);
            for f = {'KK', 'KL', 'LK', 'LL'}
                ref = m.evl.(f{1});
                tc.verifyEqual(e.(f{1}), ref, 'AbsTol', tol.(f{1})*max(abs(ref(:))), f{1});
            end
        end

        function cachedDisplacementKernelsAreReproducible(tc)
            m = tc.Ex;
            [~, d] = evalc('computeAllDisplacementKernelsBem(m.obs, m.rcv, m.shz, m.bnd, 1)');
            for f = {'KO', 'LO'}
                ref = m.devl.(f{1});
                tc.verifyEqual(d.(f{1}), ref, 'AbsTol', 1e-8*max(abs(ref(:))), f{1});
            end
        end
    end
end
