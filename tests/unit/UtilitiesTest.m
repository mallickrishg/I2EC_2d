classdef (SharedTestFixtures = {I2ECPathFixture}) UtilitiesTest < matlab.unittest.TestCase
    %UTILITIESTEST  Plotting and small helpers: plotpatch2d, plotshz2d, cline,
    % bluewhitered, drchrnd, filenameGenerator,
    % create_megathrust.  Figures are created invisibly.

    properties
        Ex
    end

    methods (TestClassSetup)
        function loadExample(tc)
            tc.Ex = loadExampleModel(false);
        end
    end

    methods (TestMethodTeardown)
        function closeFigures(~)
            close all force
        end
    end

    methods (Test)
        function plotPatchDrawsEveryPatch(tc)
            figure;
            plotpatch2d(tc.Ex.rcv);
            tc.verifyNumElements(findobj(gca, 'Type', 'line'), tc.Ex.rcv.N);
        end

        function plotPatchWithColoursUsesCline(tc)
            figure;
            plotpatch2d(tc.Ex.rcv, tc.Ex.rcv.xc(:,2)/1e3);
            tc.verifyNumElements(findobj(gca, 'Type', 'patch'), tc.Ex.rcv.N);
        end

        function plotShearZoneBothModes(tc)
            figure;
            plotshz2d(tc.Ex.shz, tc.Ex.shz.xc(:,2));
            tc.verifyNumElements(findobj(gca, 'Type', 'patch'), 1);
            figure;
            plotshz2d(tc.Ex.shz);
            tc.verifyNumElements(findobj(gca, 'Type', 'line'), tc.Ex.shz.N);
        end

        function blueWhiteRedIsAValidColormap(tc)
            figure; imagesc(peaks(20)); clim([-5 8]);
            cmap = bluewhitered(64);
            tc.verifySize(cmap, [64 3]);
            tc.verifyGreaterThanOrEqual(min(cmap(:)), 0);
            tc.verifyLessThanOrEqual(max(cmap(:)), 1);
        end

        function dirichletSamplesSumToOne(tc)
            rng(3);
            r = drchrnd([1 2 3], 50);
            tc.verifySize(r, [50 3]);
            tc.verifyEqual(sum(r, 2), ones(50, 1), 'AbsTol', 1e-12);
            tc.verifyGreaterThanOrEqual(min(r(:)), 0);
        end

        function filenameGeneratorCreatesNumberedFolders(tc)
            folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            f1 = filenameGenerator(folder);
            f2 = filenameGenerator(folder);
            tc.verifyEqual(f1, [folder '/000001/results.mat']);
            tc.verifyEqual(f2, [folder '/000002/results.mat']);
            tc.verifyTrue(isfolder(fileparts(f2)));
        end

        function createMegathrustBuildsRequestedPatches(tc)
            folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            tc.applyFixture(matlab.unittest.fixtures.CurrentFolderFixture(folder));
            dip = atand(50/200);
            rcv = create_megathrust(tc.Ex.em, 0, 10, dip, hypot(200e3, 50e3), 50);
            tc.verifyEqual(sum(rcv.W), hypot(200e3, 50e3), 'RelTol', 1e-10);
            tc.verifyEqual(rcv.dip, dip*ones(rcv.N,1), 'AbsTol', 1e-10);
            % patch count: see GeometryTest/uniformSegmentGivesUniformPatches
        end
    end

    methods (Test, TestTags = {'KnownIssue'})
        function createMegathrustLeavesWorkingFolderClean(tc)
            % writes megathrust2d.seg into the current folder: run from the
            % repository root it creates a stray file there
            folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            tc.applyFixture(matlab.unittest.fixtures.CurrentFolderFixture(folder));
            create_megathrust(tc.Ex.em, 0, 10, 14, 1e5, 10);
            tc.verifyEmpty(dir(fullfile(folder, '*.seg')), 'megathrust2d.seg written to pwd');
        end
    end
end
