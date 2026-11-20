classdef (SharedTestFixtures = {I2ECPathFixture}, TestTags = {'Slow'}) ExampleScriptsTest < matlab.unittest.TestCase
    %EXAMPLESCRIPTSTEST  End-to-end runs of the examples in examples/ (slow).
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
        function analyticNewtonianExampleRuns(tc)
            tc.verifyEmpty(ExampleScriptsTest.runScript(ExampleScriptsTest.example('example2_analyticNewtonian')));
        end

        function analyticVsNumericalExampleAgrees(tc)
            % example 4 prints the largest difference between the solvers;
            % it must be at the ode45 tolerance
            [msg, out] = ExampleScriptsTest.runScript(ExampleScriptsTest.example('example4_analyticVsNumerical'));
            tc.verifyEmpty(msg);
            d = sscanf(regexp(out, 'slip rate [^,]*', 'match', 'once'), 'slip rate %f');
            e = sscanf(regexp(out, 'strain rate [0-9.e+-]*', 'match', 'once'), 'strain rate %f');
            tc.verifyLessThan([d e], 1e-6, out);
        end

        function numericalPowerLawExampleRuns(tc)
            % five n = 3 cycles with ode45 (several minutes)
            tc.verifyEmpty(ExampleScriptsTest.runScript(ExampleScriptsTest.example('example1_numericalPowerLaw')));
        end

        function meshConvergenceExampleRuns(tc)
            % needs the cached kernels of the three meshes (an hour to build)
            cache = fullfile(i2ecRoot(), 'kernels', 'meshconvergence');
            tc.assumeTrue(all(isfile(fullfile(cache, "evl_shearzone_h" + ["40" "20" "10"] + ".h5"))), ...
                'run examples/example3_meshConvergence once to build its kernel cache');
            tc.verifyEmpty(ExampleScriptsTest.runScript(ExampleScriptsTest.example('example3_meshConvergence')));
        end

        function longTermVelocityPlotRuns(tc)
            tc.verifyEmpty(ExampleScriptsTest.runScript(ExampleScriptsTest.example('plotLongtermVelocity')));
        end

        function negateExactSolverVerificationPasses(tc)
            tc.verifyEmpty(ExampleScriptsTest.runScript(fullfile('tests', 'test_negate_exact_solver.m')));
        end
    end

    methods (Static)
        function [msg, out] = runScript(file)
            % run a script or no-argument function; return '' or the error
            % text, and the printed output
            msg = ''; out = '';
            try
                out = evalc('ExampleScriptsTest.runInOwnWorkspace(file)');
            catch err
                msg = err.getReport('basic');
            end
        end

        function file = example(name)
            file = fullfile(i2ecRoot(), 'examples', [name '.m']);
        end

        function runInOwnWorkspace(file)
            % scripts start with 'clear'; this workspace is the one cleared
            run(file);
        end
    end
end
