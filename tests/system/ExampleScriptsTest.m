classdef (SharedTestFixtures = {I2ECPathFixture}, TestTags = {'Slow'}) ExampleScriptsTest < matlab.unittest.TestCase
    %EXAMPLESCRIPTSTEST  End-to-end runs of the examples in examples/ (slow).
    %
    % Each script runs from the repository root (they use relative paths
    % such as 'examples/inputs/...') in its own function workspace, so their
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

        function transientRheologyExampleRuns(tc)
            % example 5 prints the cycle-mean departure from the long-term
            % rates for each rheology; it must be at round-off
            [msg, out] = ExampleScriptsTest.runScript(ExampleScriptsTest.example('example5_transientRheology'));
            tc.verifyEmpty(msg);
            d = str2double(regexp(out, '(?<=slip |strain )[0-9.e+-]+', 'match'));
            tc.verifyNumElements(d, 6, out);
            tc.verifyLessThan(d, 1e-6, out);
        end

        function adjointInversionExampleRecoversModel(tc)
            % example 6 checks the adjoint gradient against the tangent
            % Jacobian and recovers the KL coefficients at the noise level
            [msg, out] = ExampleScriptsTest.runScript(ExampleScriptsTest.example('example6_adjointInversion'));
            tc.verifyEmpty(msg);
            d = str2double(regexp(out, '(?<=relative difference )[0-9.e+-]+', 'match', 'once'));
            tc.verifyLessThan(d, 1e-10, out);
            chi = str2double(string(regexp(out, 'chi\^2/N = +([0-9.]+)', 'tokens')));
            err = str2double(string(regexp(out, 'xi_true\| = ([0-9.]+)', 'tokens')));
            tc.verifyLessThan(chi(end), 1.2, out);
            tc.verifyLessThan(err(end), 0.1, out);
        end

        function numericalPowerLawExampleRuns(tc)
            % five n = 3 cycles with ode45 (several minutes)
            tc.verifyEmpty(ExampleScriptsTest.runScript(ExampleScriptsTest.example('example1_numericalPowerLaw')));
        end

        function meshConvergenceExampleRuns(tc)
            % needs the cached kernels of the three meshes (about a minute to build)
            cache = fullfile(i2ecRoot(), 'examples', 'kernels', 'meshconvergence');
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
