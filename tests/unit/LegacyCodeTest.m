classdef (SharedTestFixtures = {I2ECPathFixture}) LegacyCodeTest < matlab.unittest.TestCase
    %LEGACYCODETEST  Duplicated, unused and broken files (see tasks.md §1-§3).
    %
    % Equivalence tests show whether an old copy can be deleted without
    % changing results.  Static checks scan every function file in
    % functions/, +geometry/ and meshing/.

    properties
        Ex
    end

    methods (TestClassSetup)
        function loadExample(tc)
            tc.Ex = loadExampleModel(false);
        end
    end

    methods (Test)
        %% ---------------- duplicates: same results as the kept version ----------------
        function functionsLDdispHSEqualsPackageCopy(tc)
            x = linspace(-30, 30, 7)'; y = -linspace(1, 40, 7)';
            u1 = LDdispHS(x, y, 3, -10, 2, -0.4, 1, 0.2, 0.25);
            u2 = geometry.LDdispHS(x, y, 3, -10, 2, -0.4, 1, 0.2, 0.25);
            tc.verifyEqual(u1, u2);
        end

        function functionsEdgeStressEqualsPackageCopy(tc)
            m = [5 2 10 30 1]; x = linspace(-30, 30, 7)'; z = -linspace(1, 40, 7)';
            [a1, b1, c1] = EdgeStress(m, x, z, 0.25, 30e3);
            [a2, b2, c2] = geometry.EdgeStress(m, x, z, 0.25, 30e3);
            tc.verifyEqual([a1 b1 c1], [a2 b2 c2]);
        end

        function oldDisplacementGFsEqualCurrentKernels(tc)
            rcv = tc.Ex.rcv;
            obs = [linspace(-100e3, 400e3, 9)' zeros(9,1)];
            [Gdx, Gdz] = geometry.computeFaultDisplacementKernels(rcv, obs);
            for Gd = {computedisplacementGFs(rcv, obs), geometry.computedisplacementGFs(rcv, obs)}
                tc.verifyEqual(Gd{1}(1:2:end,:), Gdx, 'AbsTol', 1e-14);
                tc.verifyEqual(Gd{1}(2:2:end,:), Gdz, 'AbsTol', 1e-14);
            end
        end

        function oldEdgeDislocationTractionKernelsEqualCurrentKernels(tc)
            % computetractionkernels (edge dislocations, 2020) vs
            % geometry.computeFaultTractionKernels (line dislocations, 2023)
            rcv = tc.Ex.rcv;
            [Kdd, Kdn] = geometry.computeFaultTractionKernels(rcv, rcv);
            for c = {@computetractionkernels, @geometry.computetractionkernels}
                [Kt, Ks] = c{1}(rcv, rcv);
                tc.verifyEqual(Kt, Kdd, 'RelTol', 1e-6, 'AbsTol', 1e-9*max(abs(Kdd(:))), func2str(c{1}));
                tc.verifyEqual(Ks, Kdn, 'RelTol', 1e-6, 'AbsTol', 1e-9*max(abs(Kdn(:))), func2str(c{1}));
            end
        end
    end

    methods (Test, TestTags = {'KnownIssue'})
        function stressFaultBemRuns(tc)
            % calls computeFullTractionKernels, which exists nowhere
            rcv = tc.Ex.rcv;
            f = @() computeStressFault_bem([-1 0 1e3 14 1], tc.Ex.bnd, rcv.xc(:,1), rcv.xc(:,2), 0.25, 30e3);
            tc.verifyWarningFree(f);
        end

        function functionNamesMatchFileNames(tc)
            % a mismatch means the function is only reachable by file name
            % (getStrainratesLongterm_old.m declares getStrainratesLongterm)
            bad = {};
            for f = LegacyCodeTest.functionFiles()'
                declared = LegacyCodeTest.declaredName(f{1});
                [~, name] = fileparts(f{1});
                if ~isempty(declared) && ~strcmp(declared, name)
                    bad{end+1} = sprintf('%s declares %s', LegacyCodeTest.relPath(f{1}), declared); %#ok<AGROW>
                end
            end
            tc.verifyEmpty(bad, strjoin(bad, newline));
        end

        function noDuplicateFileNames(tc)
            % the same function name in two folders (functions/ and +geometry/)
            files = LegacyCodeTest.functionFiles();
            names = cellfun(@(f) LegacyCodeTest.baseName(f), files, 'UniformOutput', false);
            [u, ~, j] = unique(names);
            dup = u(accumarray(j, 1) > 1);
            tc.verifyEmpty(dup, ['duplicated: ' strjoin(dup, ', ')]);
        end
    end

    methods (Static)
        function files = functionFiles()
            root = i2ecRoot();
            files = {};
            for d = {'functions', '+geometry', 'meshing'}
                l = dir(fullfile(root, d{1}, '*.m'));
                files = [files; fullfile({l.folder}, {l.name})']; %#ok<AGROW>
            end
        end

        function name = declaredName(file)
            % name after 'function [...] =' or 'classdef', '' for scripts
            name = '';
            lines = readlines(file);
            for k = 1:numel(lines)
                s = strtrim(lines(k));
                if s == "" || startsWith(s, "%")
                    continue
                end
                tok = regexp(s, '^function\s+(?:\[?[^=]*\]?\s*=\s*)?([A-Za-z]\w*)', 'tokens', 'once');
                if isempty(tok)
                    tok = regexp(s, '^classdef\s*(?:\([^)]*\))?\s*([A-Za-z]\w*)', 'tokens', 'once');
                end
                if ~isempty(tok)
                    name = char(tok{1});
                end
                return
            end
        end

        function n = baseName(f)
            [~, n] = fileparts(f);
        end

        function p = relPath(f)
            p = erase(f, [i2ecRoot() filesep]);
        end
    end
end
