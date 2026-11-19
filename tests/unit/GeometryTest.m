classdef (SharedTestFixtures = {I2ECPathFixture}) GeometryTest < matlab.unittest.TestCase
    %GEOMETRYTEST  Mesh classes in +geometry and the half-space kernel wrappers.

    properties
        Ex   % example model (inputs/), geometry only
    end

    methods (TestClassSetup)
        function loadExample(tc)
            tc.Ex = loadExampleModel(false);
        end
    end

    methods (Test)
        %% ---------------- earth model ----------------
        function earthModelStoresParameters(tc)
            em = geometry.LDhs(30e3, 0.25);
            tc.verifyEqual([em.G em.nu], [30e3 0.25]);
            tc.verifyTrue(isa(em, 'geometry.earthModel'));
        end

        function earthModelRejectsBadPoissonRatio(tc)
            tc.verifyError(@() geometry.LDhs(30e3, 0.6), ?MException);
            tc.verifyError(@() geometry.LDhs(30e3, -1.5), ?MException);
            tc.verifyError(@() geometry.LDhs(-1, 0.25), ?MException);
        end

        %% ---------------- fault receiver ----------------
        function megathrustPatchCountAndLength(tc)
            rcv = tc.Ex.rcv;
            tc.verifyGreaterThanOrEqual(rcv.N, 49);
            tc.verifyEqual(sum(rcv.W), 206155.281280883, 'RelTol', 1e-10, 'patches cover the segment');
            tc.verifyEqual(rcv.dip, 14.036243467926479*ones(rcv.N,1), 'AbsTol', 1e-10);
        end

        function receiverUnitVectorsAreOrthonormal(tc)
            for r = {tc.Ex.rcv, tc.Ex.bnd}
                o = r{1};
                tc.verifyEqual(vecnorm(o.dv, 2, 2), ones(o.N,1), 'AbsTol', 1e-12);
                tc.verifyEqual(vecnorm(o.nv, 2, 2), ones(o.N,1), 'AbsTol', 1e-12);
                tc.verifyEqual(sum(o.dv .* o.nv, 2), zeros(o.N,1), 'AbsTol', 1e-12);
            end
        end

        function receiverPatchesAreContiguous(tc)
            % the bottom end of patch k is the top end of patch k+1
            rcv = tc.Ex.rcv;
            bottom = rcv.x - rcv.W .* rcv.dv;
            tc.verifyEqual(rcv.x(2:end,:), bottom(1:end-1,:), 'AbsTol', 1e-6);
            tc.verifyEqual(rcv.xc, rcv.x - rcv.W/2 .* rcv.dv, 'AbsTol', 1e-9);
        end

        function receiverDepthsAreBelowSurface(tc)
            tc.verifyLessThanOrEqual(max(tc.Ex.rcv.xc(:,2)), 0);
            tc.verifyLessThanOrEqual(max(tc.Ex.shz.xc(:,2)), 0);
        end

        function boundaryIsDirichletWithZeroDisplacement(tc)
            b = tc.Ex.bnd;
            tc.verifyEqual(b.N, 112);
            tc.verifyTrue(all(b.Vpl == 0), 'boundary2d.seg: every patch u = 0');
            tc.verifyEqual([b.Vx; b.Vz], zeros(2*b.N, 1));
        end

        function receiverRejectsMissingFile(tc)
            tc.verifyError(@() geometry.receiver('no_such_file.seg', tc.Ex.em), ?MException);
        end

        %% ---------------- shear zone receiver ----------------
        function shearZoneMeshIsValid(tc)
            shz = tc.Ex.shz;
            tc.verifyEqual(shz.N, 521);
            tc.verifyGreaterThan(min(shz.area), 0, 'no degenerate triangles');
            tc.verifyEqual(shz.xc, (shz.A + shz.B + shz.C)/3, 'AbsTol', 1e-6);
            tc.verifyEqual(size(shz.tri), [shz.N 3]);
            tc.verifyEqual(max(shz.tri(:)), size(shz.vert, 1), 'all vertices referenced');
        end

        function shearZoneAreaEqualsMeshArea(tc)
            m = makeSmallModel(tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder);
            tc.verifyEqual(sum(m.shz.area), 100e3 * 60e3, 'RelTol', 1e-12);
        end

        %% ---------------- half-space kernel wrappers ----------------
        function faultTractionKernelIsReciprocal(tc)
            % W_i K(i,j) = W_j K(j,i) for well-separated patches (Betti)
            rcv = tc.Ex.rcv;
            [Kdd, ~, ~, Knn] = geometry.computeFaultTractionKernels(rcv, rcv);
            d = sqrt((rcv.xc(:,1) - rcv.xc(:,1)').^2 + (rcv.xc(:,2) - rcv.xc(:,2)').^2);
            far = d > 3*max(rcv.W);
            for K = {Kdd, Knn}
                M = K{1} .* rcv.W;
                r = reciprocityRatio(M, M.', far);
                tc.verifyEqual(median(r), 1, 'AbsTol', 1e-3);
                % midpoint collocation: ~2% spread at 3 patch widths
                tc.verifyEqual(prctile(r, [5 95]), [1 1], 'AbsTol', 5e-2);
            end
            tc.verifyTrue(all(diag(Kdd) < 0), 'self shear stress opposes slip');
        end

        function faultStressKernelsAgreeWithTractionKernels(tc)
            % computeFaultStressKernels returns -Sxz (depth-positive frame);
            % rotating it back must reproduce computeFaultTractionKernels.
            rcv = tc.Ex.rcv;
            [Kdd, Kdn] = geometry.computeFaultTractionKernels(rcv, rcv);
            [Sxx, Szz, mSxz] = geometry.computeFaultStressKernels(rcv, rcv);
            Sxz = -mSxz;
            tx = Sxx .* rcv.nv(:,1) + Sxz .* rcv.nv(:,2);
            tz = Sxz .* rcv.nv(:,1) + Szz .* rcv.nv(:,2);
            tc.verifyEqual(rcv.dv(:,1).*tx + rcv.dv(:,2).*tz, Kdd, 'RelTol', 1e-10, 'AbsTol', 1e-12);
            tc.verifyEqual(rcv.nv(:,1).*tx + rcv.nv(:,2).*tz, Kdn, 'RelTol', 1e-10, 'AbsTol', 1e-12);
        end

        function faultDisplacementKernelJumpsByUnitSlip(tc)
            % straddle patch 20 at its centre: uniform unit slip on the
            % whole fault gives a displacement jump of 1 along the fault
            rcv = tc.Ex.rcv;
            k = 20;
            d = 1;   % m
            obs = [rcv.xc(k,:) + d*rcv.nv(k,:); rcv.xc(k,:) - d*rcv.nv(k,:)];
            [Gdx, Gdz] = geometry.computeFaultDisplacementKernels(rcv, obs);
            jump = [sum(Gdx(1,:) - Gdx(2,:)), sum(Gdz(1,:) - Gdz(2,:))];
            tc.verifyEqual(norm(jump), 1, 'AbsTol', 1e-3);
            tc.verifyEqual(abs(jump * rcv.dv(k,:)'), 1, 'AbsTol', 1e-3, 'jump is along the dip direction');
        end
    end

    methods (Test, TestTags = {'KnownIssue'})
        function faultDisplacementKernelOnPatchCentre(tc)
            % An observation point exactly on a patch centre is nudged along
            % rcv.nv(r==0,:), but r==0 is a mask over observation points,
            % not patches.  With Nobs ~= N this takes another patch's normal,
            % or errors when the point's row is beyond rcv.N; the nudge
            % should use rcv.nv(k,:).  Here the coincident point is row 61.
            rcv = tc.Ex.rcv;
            obs = [linspace(-100e3, 400e3, 60)' zeros(60,1); rcv.xc(5,:)];
            f = @() geometry.computeFaultDisplacementKernels(rcv, obs);
            tc.verifyWarningFree(f);
            [Gdx, Gdz] = geometry.computeFaultDisplacementKernels(rcv, obs);
            tc.verifyTrue(all(isfinite([Gdx(:); Gdz(:)])));
        end

        function uniformSegmentGivesUniformPatches(tc)
            % megathrust2d.seg asks for W0 = W/50 with growth factor 1.
            % flt2flt merges the last two patches when the remainder is
            % exactly 2*W0 (floating-point '>' on equal values), so the
            % mesh has 49 patches and the deepest one is twice as wide.
            % create_megathrust(..., Nfault) returns Nfault-1 patches.
            W = tc.Ex.rcv.W;
            tc.verifyEqual(numel(W), 50, 'number of patches');
            tc.verifyEqual(max(W) / min(W), 1, 'RelTol', 1e-9, 'width ratio deepest/other');
        end
    end
end
