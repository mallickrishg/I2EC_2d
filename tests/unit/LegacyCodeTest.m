classdef (SharedTestFixtures = {I2ECPathFixture}) LegacyCodeTest < matlab.unittest.TestCase
    %LEGACYCODETEST  Static checks on the code base (see tasks.md §1-§3).
    %
    % Scans every function file in functions/ and +geometry/.

    methods (Test)
        function functionNamesMatchFileNames(tc)
            % a mismatch means the function is only reachable by file name
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
            for d = {'functions', '+geometry'}
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
