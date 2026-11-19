classdef (SharedTestFixtures = {I2ECPathFixture}, TestTags = {'Slow'}) ExampleScriptsTest < matlab.unittest.TestCase
    %EXAMPLESCRIPTSTEST  End-to-end runs of the top-level scripts (slow).
    %
    % Each script runs from the repository root (they use relative paths
    % such as 'inputs/...') in its own function workspace, so their
    % 'clear' does not touch the test.  Figures stay invisible.

    methods (TestMethodSetup)
        function goToRoot(tc)
            tc.applyFixture(matlab.unittest.fixtures.CurrentFolderFixture(i2ecRoot()));
        end
    end

    methods (TestMethodTeardown)
        function closeFigures(~)
            close all force
        end
    end

    methods (Test)
        function analyticSingleCycleScriptRuns(tc)
            tc.verifyEmpty(ExampleScriptsTest.runScript('MAIN_analyticsolve_singleeqcycle.m'));
        end

        function longTermVelocityPlotRuns(tc)
            tc.verifyEmpty(ExampleScriptsTest.runScript('PLOT_longtermVelocity.m'));
        end

        function negateExactSolverVerificationPasses(tc)
            tc.verifyEmpty(ExampleScriptsTest.runScript(fullfile('tests', 'test_negate_exact_solver.m')));
        end

        function numericalSingleCycleScriptRuns(tc)
            % rebuilds all kernels on the full mesh (~2 min)
            tc.verifyEmpty(ExampleScriptsTest.runScript('MAIN_i2ec_singleeqcycle.m'));
        end
    end

    methods (Static)
        function msg = runScript(file)
            % run a script or no-argument function; return '' or the error text
            msg = '';
            try
                evalc('ExampleScriptsTest.runInOwnWorkspace(file)');
            catch err
                msg = err.getReport('basic');
            end
        end

        function runInOwnWorkspace(file)
            % scripts start with 'clear'; this workspace is the one cleared
            run(file);
        end
    end
end
