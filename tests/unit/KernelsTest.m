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

    properties (TestParameter)
        n = {1, 3}
        mode = struct('zero', 1, 'damp', 2)
        quadrature = {3, 7}
    end

    properties
        Sm     % small model struct (makeSmallModel)
        Evl    % small-model centroid kernels ('cellaverage', 1, not certified, not modified)
        EvlCert % small-model default kernels (3-point cell average, certified)
        Devl   % small-model displacement kernels
        Ex     % example model with cached evl_orig / devl
    end

    methods (TestClassSetup)
        function buildKernels(tc)
            folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            tc.Sm = makeSmallModel(folder);
            [~, tc.Evl] = evalc('computeAllStressKernelsBem(tc.Sm.rcv, tc.Sm.shz, tc.Sm.bnd, ''cellaverage'', 1, ''kernelcertify'', false)');
            [~, tc.EvlCert] = evalc('computeAllStressKernelsBem(tc.Sm.rcv, tc.Sm.shz, tc.Sm.bnd)');
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
            [~, evm] = evalc('computeAllStressKernelsBem(tc.Sm.rcv, tc.Sm.shz, tc.Sm.bnd, ''cellaverage'', 1, ''kernelcertify'', false, ''kernelmodify'', 1)');
            L = [evm.LL(:,:,1,1) evm.LL(:,:,1,2); evm.LL(:,:,2,1) evm.LL(:,:,2,2)];
            scale = max(abs(eig(L)));
            tc.verifyLessThanOrEqual(max(real(eig(L))) / scale, 1e-10);
            tc.verifyEqual(evm.KK, tc.Evl.KK, 'other blocks unchanged');
            tc.verifyEqual(evm.LK, tc.Evl.LK);
        end

        function dampedKernelHasOnlyDecayingModes(tc)
            % kernelmodify = 2 on a centroid kernel: each growing
            % eigenvalue becomes -max|lambda|
            e0 = tc.Evl;
            ed = modifyGrowingModes(e0, 2);
            lam0 = eig(KernelsTest.deviatoricBlock(e0));
            lam = eig(KernelsTest.deviatoricBlock(ed));
            scale = max(abs(lam0));
            tc.verifyLessThan(max(real(lam)) / scale, 0, 'no growing or neutral mode');
            nGrowing = sum(real(lam0) > 0);
            atMax = @(x) abs(x + scale) < 1e-6*scale;   % the kernel's own -max|lambda| counts too
            tc.verifyEqual(sum(atMax(lam)), nGrowing + sum(atMax(lam0)), ...
                'each growing mode now decays at -max|lambda|');
            tc.verifyEqual(sort(real(lam(~atMax(lam)))), ...
                sort(real(lam0(real(lam0) <= 0 & ~atMax(lam0)))), 'AbsTol', 1e-8*scale, 'other modes unchanged');
            for f = {'KK', 'KL', 'LK'}
                tc.verifyEqual(ed.(f{1}), e0.(f{1}), f{1});
            end
        end

        function modifiedKernelsHaveNoGrowingModesForPowerLaw(tc, n, mode)
            % Linearise the mantle system ds/dt = L * edot(s) about random
            % stress states: d(edot)/ds = alpha tau^(n-1) (I + (n-1) s s'/tau^2)
            % is positive definite per cell, so the number of growing
            % modes is that of L for any n; a modified kernel has none.
            L = KernelsTest.deviatoricBlock(modifyGrowingModes(tc.Evl, mode));
            M = tc.Sm.shz.N;
            rng(3);
            for trial = 1:2
                tau = 10.^(-1 + 2*rand(M, 1)); ang = 2*pi*rand(M, 1);
                s22 = tau.*cos(ang); s23 = tau.*sin(ang);
                c = 1e-14*(1 + 9*rand(M, 1)) .* tau.^(n-1);
                D11 = c.*(1 + (n-1)*s22.^2./tau.^2);
                D22 = c.*(1 + (n-1)*s23.^2./tau.^2);
                D12 = c.*(n-1).*s22.*s23./tau.^2;
                lam = eig(L * [diag(D11) diag(D12); diag(D12) diag(D22)]);
                tc.verifyLessThan(max(real(lam)) / max(abs(lam)), 1e-10, sprintf('trial %d', trial));
            end
        end

        function kernelModifyOptionMatchesModifyGrowingModes(tc)
            [~, ed] = evalc('computeAllStressKernelsBem(tc.Sm.rcv, tc.Sm.shz, tc.Sm.bnd, ''cellaverage'', 1, ''kernelcertify'', false, ''kernelmodify'', 2)');
            ref = modifyGrowingModes(tc.Evl, 2);
            tc.verifyEqual(ed.LL, ref.LL, 'AbsTol', 1e-10*max(abs(ref.LL(:))));
            tc.verifyEqual(ed.KK, tc.Evl.KK);
        end

        function cellAveragedMantleKernelHasNoGrowingModes(tc, quadrature)
            % Averaging the stress over each receiving cell (Galerkin)
            % instead of sampling it at the centroid removes the growing
            % modes at the source and restores reciprocity of LL.
            m = tc.Sm;
            [~, e] = evalc('computeAllStressKernelsBem(m.rcv, m.shz, m.bnd, ''kernelmodify'', 0, ''cellaverage'', quadrature)');
            L = KernelsTest.deviatoricBlock(e);
            lam = eig(L);
            tc.verifyEqual(sum(real(lam) > 0), 0, sprintf('max Re/max|lambda| = %.2e', max(real(lam))/max(abs(lam))));
            S = diag([m.shz.area; m.shz.area]) * L;
            tc.verifyLessThan(norm(S - S', 'fro') / norm(S, 'fro'), 1e-12, 'certified: area-weighted LL is symmetric');
            tc.verifyLessThanOrEqual(max(eig((S + S')/2)), 1e-12*max(abs(eig(S))), 'certified: negative semidefinite');
            tc.verifyLessThan(e.certificate.asymmetry, 0.05, 'quadrature asymmetry removed by certification');
            tc.verifyLessThanOrEqual(e.certificate.maxEigenvalue, e.certificate.tolerance);
            tc.verifyEqual(m.shz.xc, (m.shz.A + m.shz.B + m.shz.C)/3, 'AbsTol', 1e-6, 'centroids restored');
            tc.verifyEqual(e.KK, tc.Evl.KK, 'fault-fault block unchanged');
            tc.verifyEqual(e.LK, tc.Evl.LK, 'mantle-to-fault block unchanged');
        end

        function certificationRejectsCentroidKernel(tc)
            % centroid kernels are not negative semidefinite (the
            % example mesh: 65 positive eigenvalues in the symmetric part)
            tc.verifyError(@() certifyMantleKernel(tc.Evl, tc.Sm.shz), 'certifyMantleKernel:notNegative');
        end

        function certifiedKernelHasNoGrowingModesForPowerLaw(tc, n)
            % the guarantee: no growing modes for any SPD power-law
            % Jacobian (random stress states, magnitudes and prefactors)
            m = tc.Sm;
            [~, e] = evalc('computeAllStressKernelsBem(m.rcv, m.shz, m.bnd, ''kernelmodify'', 0, ''cellaverage'', 3)');
            L = KernelsTest.deviatoricBlock(e);
            M = m.shz.N;
            rng(5);
            for trial = 1:5
                tau = 10.^(-2 + 4*rand(M, 1)); ang = 2*pi*rand(M, 1);
                s22 = tau.*cos(ang); s23 = tau.*sin(ang);
                c = 10.^(-16 + 4*rand(M, 1)) .* tau.^(n-1);
                D11 = c.*(1 + (n-1)*s22.^2./tau.^2);
                D22 = c.*(1 + (n-1)*s23.^2./tau.^2);
                D12 = c.*(n-1).*s22.*s23./tau.^2;
                lam = eig(L * [diag(D11) diag(D12); diag(D12) diag(D22)]);
                tc.verifyLessThanOrEqual(max(real(lam)), 1e-10*max(abs(lam)), sprintf('trial %d', trial));
            end
        end

        function centroidCollocationHasGrowingModesOnTestMesh(tc)
            % the reference for the test above: the default (centroid)
            % kernel of the same mesh does have growing modes
            tc.verifyGreaterThan(sum(real(eig(KernelsTest.deviatoricBlock(tc.Evl))) > 0), 0);
        end

        function defaultIsCertifiedThreePointAverage(tc)
            m = tc.Sm;
            [~, e] = evalc('computeAllStressKernelsBem(m.rcv, m.shz, m.bnd, ''cellaverage'', 3, ''kernelcertify'', true, ''kernelmodify'', 0)');
            tc.verifyEqual(tc.EvlCert, e);
            tc.verifyTrue(isfield(tc.EvlCert, 'certificate'));
            tc.verifyFalse(isfield(tc.Evl, 'certificate'), 'centroid kernels are not certified');
        end

        function cachedExampleKernelIsCertified(tc)
            % kernels/evl_orig.h5 holds the default (certified) kernels:
            % no growing modes, area-weighted LL exactly symmetric
            e = tc.Ex.evl;
            tc.verifyTrue(isfield(e, 'certificate'), 'kernels/evl_orig.h5 is not a certified kernel');
            L = KernelsTest.deviatoricBlock(e);
            lam = eig(L);
            tc.verifyLessThanOrEqual(max(real(lam)), 1e-12*max(abs(lam)), ...
                sprintf('%d eigenvalues with Re > 0', sum(real(lam) > 0)));
            S = [tc.Ex.shz.area; tc.Ex.shz.area] .* L;
            tc.verifyLessThan(norm(S - S', 'fro')/norm(S, 'fro'), 1e-12);
        end

        function kernelFilesRoundTrip(tc)
            % saveKernels / loadKernels: exact round trip of the stress and
            % displacement kernels, nested structs and geometry
            folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            m = tc.Sm;
            f1 = fullfile(folder, 'evl.h5'); f2 = fullfile(folder, 'devl.h5');
            saveKernels(f1, tc.EvlCert, 'rcv', m.rcv, 'shz', m.shz, 'description', 'test');
            saveKernels(f2, tc.Devl, 'obs', m.obs);
            [e, g, info] = loadKernels(f1);
            tc.verifyEqual(e, tc.EvlCert, 'stress kernels, certificate and options');
            tc.verifyEqual(loadKernels(f2), tc.Devl, 'displacement kernels');
            tc.verifyEqual(g.mantle.vertices, m.shz.vert);
            tc.verifyEqual(double(g.mantle.triangles), m.shz.tri);
            tc.verifyEqual(g.fault.W, m.rcv.W);
            tc.verifyEqual(info.format, 'I2EC kernels');
            tc.verifyEqual(info.description, 'test');
            [~, g2] = loadKernels(f2);
            tc.verifyEqual(g2.obs, m.obs);
        end

        function kernelFilesShowTheMatlabShapeToRowMajorReaders(tc)
            % h5info reports dimensions fastest-first; a row-major reader
            % (h5py/NumPy) sees them reversed, which must be the MATLAB
            % size, so f['KL'][i, j, c] = KL(i+1, j+1, c+1)
            folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            f = fullfile(folder, 'evl.h5');
            saveKernels(f, tc.EvlCert);
            for name = {'KK', 'KL', 'LK', 'LL'}
                I = h5info(f, ['/' name{1}]);
                tc.verifyEqual(fliplr(I.Dataspace.Size), size(tc.EvlCert.(name{1})), name{1});
                tc.verifyFalse(isempty(h5readatt(f, ['/' name{1}], 'description')));
            end
            % one element read through the low-level, C-order interface
            fid = H5F.open(f); did = H5D.open(fid, '/KL');
            sid = H5D.get_space(did);
            H5S.select_hyperslab(sid, 'H5S_SELECT_SET', [2 4 1], [], [1 1 1], []);   % C order, 0-based
            msid = H5S.create_simple(1, 1, []);
            v = H5D.read(did, 'H5ML_DEFAULT', msid, sid, 'H5P_DEFAULT');
            H5S.close(msid); H5S.close(sid); H5D.close(did); H5F.close(fid);
            tc.verifyEqual(v, tc.EvlCert.KL(3, 5, 2), 'KL[2,4,1] (0-based, row-major) = KL(3,5,2)');
        end

        function loadKernelsRejectsOtherFiles(tc)
            folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            tc.verifyError(@() loadKernels(fullfile(folder, 'none.h5')), 'loadKernels:missing');
            f = fullfile(folder, 'other.h5');
            h5create(f, '/x', 3); h5write(f, '/x', [1 2 3]');
            tc.verifyError(@() loadKernels(f), 'loadKernels:format');
        end

        function invalidKernelModifyIsRejected(tc)
            m = tc.Sm;
            tc.verifyError(@() modifyGrowingModes(tc.Evl, 3), 'modifyGrowingModes:mode');
            tc.verifyError(@() computeAllStressKernelsBem(m.rcv, m.shz, m.bnd, 'kernelmodify', 3), ?MException);
            tc.verifyError(@() computeAllStressKernelsBem(m.rcv, m.shz, m.bnd, 'cellaverage', 2), ?MException);
            tc.verifyError(@() computeAllStressKernelsBem(m.rcv, m.shz, m.bnd, 'cellaverage', 1), ?MException, ...
                'certification of a centroid kernel must be refused');
            tc.verifyError(@() computeAllStressKernelsBem(m.rcv, m.shz, m.bnd, 'kernelcertify', 'yes'), ?MException);
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
            % kernels/*.h5 must be built on the meshes in examples/inputs/
            e = tc.Ex.evl;
            tc.verifySize(e.KK, [tc.Ex.rcv.N tc.Ex.rcv.N]);
            tc.verifySize(e.LL, [tc.Ex.shz.N tc.Ex.shz.N 2 2]);
            tc.verifySize(tc.Ex.devl.KO, [size(tc.Ex.obs,1) tc.Ex.rcv.N 2]);
        end
    end

    methods (Test)
        function stressKernelsIgnoreBoundaryVelocity(tc)
            % The BEM kernels solve the homogeneous problem: the prescribed
            % boundary motion (bcx/bcz in the .seg file -> boundary.Vx, Vz)
            % is not used, so every column is linear in its unit source.
            % (It used to enter the Dirichlet data of every column.)
            % own copy of the model: receiver objects are handles
            m = makeSmallModel(tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder);
            [K0, S0] = computeFaultTractionKernelsBem(m.rcv, m.rcv, m.bnd);
            L0 = computeShzStressKernelsBem(m.shz, m.shz, m.bnd);
            T0 = computeShzStressKernelsBem(m.shz, m.rcv, m.bnd);
            % non-uniform (a uniform value is a rigid translation: no stress)
            m.bnd.Vx = linspace(0, 1, m.bnd.N)';
            m.bnd.Vz = linspace(1, -1, m.bnd.N)'.^2;
            [K1, S1] = computeFaultTractionKernelsBem(m.rcv, m.rcv, m.bnd);
            L1 = computeShzStressKernelsBem(m.shz, m.shz, m.bnd);
            T1 = computeShzStressKernelsBem(m.shz, m.rcv, m.bnd);
            tc.verifyEqual(K1, K0, 'fault-fault shear traction');
            tc.verifyEqual(S1, S0, 'fault-fault normal traction');
            tc.verifyEqual(L1, L0, 'mantle-mantle stress');
            tc.verifyEqual(T1, T0, 'mantle-fault traction');
        end
    end

    methods (Test, TestTags = {'Slow'})
        function cachedStressKernelsAreReproducible(tc)
            % rebuild evl_orig (default: certified 3-point) from examples/inputs/ and
            % compare to kernels/evl_orig.h5
            % (minutes; catches a stale cache after a mesh or kernel change)
            m = tc.Ex;
            [~, e] = evalc('computeAllStressKernelsBem(m.rcv, m.shz, m.bnd)');
            % every block is closed form (line dislocations for KK, KL; the
            % differentiated triangle displacement for LK, LL, since
            % 2026-11-20), so the cache reproduces to rounding. A cache built
            % with the earlier finite-difference triangle stress differs by
            % ~1e-5 in LK and 4e-5 in LL and fails here
            tol = struct('KK', 1e-8, 'KL', 1e-8, 'LK', 1e-8, 'LL', 1e-8);
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

    methods (Static)
        function L = deviatoricBlock(e)
            L = [e.LL(:,:,1,1) e.LL(:,:,1,2); e.LL(:,:,2,1) e.LL(:,:,2,2)];
        end
    end
end
