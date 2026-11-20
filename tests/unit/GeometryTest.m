classdef (SharedTestFixtures = {I2ECPathFixture}) GeometryTest < matlab.unittest.TestCase
    %GEOMETRYTEST  Mesh classes in +geometry and the half-space kernel wrappers.

    properties
        Ex   % example model (examples/inputs/), geometry only
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

        function createdMeshIsReadableAndFillsTheOutline(tc)
            % meshing/createShearZoneMesh: PDE Toolbox mesh written as .dat
            % files that geometry.shearZoneReceiver reads; the triangles
            % tile both polygons exactly
            tc.assumeTrue(license('test', 'PDE_Toolbox') && exist('createpde', 'file') > 0, ...
                'needs the PDE Toolbox');
            folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            [vert, tri] = createShearZoneMesh('Hmax', 40, 'OutputFolder', folder, 'Name', 'mesh');
            shz = geometry.shearZoneReceiver(fullfile(folder, 'mesh'), tc.Ex.em);
            tc.verifyEqual(shz.N, size(tri, 1));
            tc.verifyEqual(shz.vert, vert, 'AbsTol', 1e-6);
            tc.verifyGreaterThan(min(shz.area), 0, 'no degenerate triangles');
            % default outline: same as the example mesh in examples/inputs/
            oc = [-100 -20; -100 -200; 450 -200; 450 -132.5; 0 -20];
            we = [201 -50; 450 -112.5; 450 -35; 201 -35];
            outline = (polyarea(oc(:,1), oc(:,2)) + polyarea(we(:,1), we(:,2))) * 1e6;
            tc.verifyEqual(sum(shz.area), outline, 'RelTol', 1e-10, 'mesh fills the polygons');
            tc.verifyEqual(sum(tc.Ex.shz.area), outline, 'RelTol', 1e-10, 'default outline = example mesh');
            % every triangle lies in one of the two polygons
            inO = inpolygon(shz.xc(:,1)/1e3, shz.xc(:,2)/1e3, oc(:,1), oc(:,2));
            inW = inpolygon(shz.xc(:,1)/1e3, shz.xc(:,2)/1e3, we(:,1), we(:,2));
            tc.verifyTrue(all(xor(inO, inW)));
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

        function faultDisplacementKernelOnPatchCentre(tc)
            % A point exactly on the centre of patch k must be nudged along
            % rcv.nv(k,:), giving the value on the +nv side.  Until
            % 2026-11-19 the normal was rcv.nv(r==0,:), indexed by the
            % observation point: here the point is row 61 > rcv.N, which
            % errored.
            rcv = tc.Ex.rcv;
            k = 5;
            obs = [linspace(-100e3, 400e3, 60)' zeros(60,1); rcv.xc(k,:)];
            f = @() geometry.computeFaultDisplacementKernels(rcv, obs);
            tc.verifyWarningFree(f);
            [Gdx, Gdz] = geometry.computeFaultDisplacementKernels(rcv, obs);
            tc.verifyTrue(all(isfinite([Gdx(:); Gdz(:)])));
            [Px, Pz] = geometry.computeFaultDisplacementKernels(rcv, rcv.xc(k,:) + 1e-6*rcv.W(k)*rcv.nv(k,:));
            tc.verifyEqual([Gdx(end,:); Gdz(end,:)], [Px; Pz], 'AbsTol', 1e-5, '+nv side of patch k');
        end

        function faultDisplacementKernelOnPatchCentreNegativeX(tc)
            % The step used to be 1e-9*mean(ox(:,1)), which is negative
            % when the observation points lie mostly at negative x, so the
            % point landed on the -nv side and the self-displacement of
            % patch k had the wrong sign.
            rcv = tc.Ex.rcv;
            k = 5;
            obs = [-300e3 0; rcv.xc(k,:)];
            [Gdx, Gdz] = geometry.computeFaultDisplacementKernels(rcv, obs);
            [Px, Pz] = geometry.computeFaultDisplacementKernels(rcv, rcv.xc(k,:) + 1e-6*rcv.W(k)*rcv.nv(k,:));
            tc.verifyEqual([Gdx(2,:); Gdz(2,:)], [Px; Pz], 'AbsTol', 1e-5, '+nv side of patch k');
        end
    end

    methods (Test)
        function uniformSegmentMergesTheLastTwoPatches(tc)
            % Known quirk, kept on purpose (won't fix, see tasks.md):
            % megathrust2d.seg asks for W0 = W/50 with growth factor 1, but
            % flt2flt merges the last two patches when the remainder is
            % 2*W0 (round-off makes 'Wt > Wc/2' true), so the example mesh
            % has 49 patches and the deepest one is 2*W0 wide. The cached
            % kernels in examples/kernels/ are built on this mesh, so this test pins
            % it: a change here means the kernels must be rebuilt.
            rcv = tc.Ex.rcv;
            W = rcv.W;
            W0 = min(W);
            tc.verifyEqual(numel(W), 49, 'number of patches');
            tc.verifyEqual(W(1:end-1), W0*ones(48,1), 'RelTol', 1e-12, 'uniform patches');
            tc.verifyEqual(W(end), 2*W0, 'RelTol', 1e-9, 'deepest patch');
            tc.verifyEqual(sum(W), 50*W0, 'RelTol', 1e-12, 'total length');
            [~, deepest] = max(rcv.xc(:,2) .* -1);
            tc.verifyEqual(deepest, numel(W), 'the wide patch is the deepest');
        end
    end
end
