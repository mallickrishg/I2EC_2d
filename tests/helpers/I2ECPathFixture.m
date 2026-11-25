classdef I2ECPathFixture < matlab.unittest.fixtures.Fixture
    %I2ECPATHFIXTURE  Put the I2EC code on the path for a test class.
    %
    % Adds, in this order (last added is searched first):
    %   <root>/functions, <root>/inverse, <root> (for the +geometry package)
    % and restores the previous path on teardown.  Also hides figures and
    % starts no parallel pool (parfor runs serially without one).
    %
    % Use with  TestClassSetup / SharedTestFixtures:
    %   classdef (SharedTestFixtures = {I2ECPathFixture}) MyTest < matlab.unittest.TestCase

    properties (Access = private)
        OldPath
        OldFigVisible
    end

    methods
        function setup(fixture)
            fixture.OldPath = path;
            root = i2ecRoot();
            addpath(fullfile(root, 'functions'));
            addpath(fullfile(root, 'inverse'));
            addpath(root);
            addpath(fullfile(root, 'tests', 'helpers'));
            fixture.OldFigVisible = get(groot, 'DefaultFigureVisible');
            set(groot, 'DefaultFigureVisible', 'off');
            fixture.SetupDescription = 'Add I2EC functions/ and +geometry to the path';
        end

        function teardown(fixture)
            path(fixture.OldPath);
            set(groot, 'DefaultFigureVisible', fixture.OldFigVisible);
            close all force
        end
    end

    methods (Access = protected)
        function tf = isCompatible(~, ~)
            tf = true;
        end
    end
end
