function varargout = runTests(mode)
%RUNTESTS  Run the I2EC test suite and summarise what works and what does not.
%
%   runTests            fast tests: tests/unit, without tag 'Slow' (~1-2 min)
%   runTests('all')     also tests/system and every 'Slow' test (rebuilds
%                       the example kernels and runs the MAIN scripts; long)
%   report = runTests(...)    struct with the test names in each group
%
% Run from anywhere; e.g. from the repository root:
%   addpath tests; runTests
% or headless:
%   matlab -batch "addpath tests; runTests"
%
% Tests tagged 'KnownIssue' check behaviour that is currently wrong
% (they are expected to fail until the code is fixed).  The summary lists:
%   - unexpected failures   (a test without 'KnownIssue' failed: regression)
%   - known issues          (still failing, as documented in tasks.md)
%   - fixed known issues    (a 'KnownIssue' test now passes: remove the tag)
%
% Writes tests/results/summary.md and tests/results/junit.xml.
% Errors at the end if there is any unexpected failure, so a CI job or
% 'matlab -batch' exits non-zero only on regressions.

import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.selectors.HasTag
import matlab.unittest.plugins.XMLPlugin

if nargin < 1
    mode = 'fast';
end
mode = validatestring(mode, {'fast', 'all'});

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'helpers'));
outdir = fullfile(here, 'results');
if ~isfolder(outdir)
    mkdir(outdir);
end

suite = TestSuite.fromFolder(fullfile(here, 'unit'));
if strcmp(mode, 'all')
    suite = [suite, TestSuite.fromFolder(fullfile(here, 'system'))];
else
    suite = suite.selectIf(~HasTag('Slow'));
end

runner = TestRunner.withTextOutput('OutputDetail', 1);
runner.addPlugin(XMLPlugin.producingJUnitFormat(fullfile(outdir, 'junit.xml')));
runner.addPlugin(matlab.unittest.plugins.DiagnosticsRecordingPlugin);
results = runner.run(suite);

% ---- classify ----
known = arrayfun(@(t) any(strcmp(t.Tags, 'KnownIssue')), suite);
failed = [results.Failed] | [results.Incomplete];
names = {results.Name};

report = struct();
report.mode = mode;
report.passed = names(~failed & ~known);
report.unexpectedFailures = names(failed & ~known);
report.knownIssues = names(failed & known);
report.fixedKnownIssues = names(~failed & known);
report.results = results;

% ---- summary ----
lines = {};
lines{end+1} = sprintf('# I2EC test summary (%s, %s)', mode, char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
lines{end+1} = '';
lines{end+1} = sprintf('- passed: %d', numel(report.passed));
lines{end+1} = sprintf('- unexpected failures: %d', numel(report.unexpectedFailures));
lines{end+1} = sprintf('- known issues still failing: %d', numel(report.knownIssues));
lines{end+1} = sprintf('- known issues now passing: %d', numel(report.fixedKnownIssues));
lines = [lines, section('Unexpected failures (regressions)', report.unexpectedFailures, results)];
lines = [lines, section('Known issues (expected to fail, see tasks.md)', report.knownIssues, results)];
lines = [lines, section('Known issues now passing (remove the KnownIssue tag)', report.fixedKnownIssues, results)];
lines{end+1} = '';
lines{end+1} = sprintf('Total time: %.1f s', sum([results.Duration]));

text = strjoin(lines, newline);
fprintf('\n%s\n', text);
writelines(string(text), fullfile(outdir, 'summary.md'));

if nargout > 0
    varargout{1} = report;
end

if ~isempty(report.unexpectedFailures)
    error('runTests:failures', '%d unexpected test failure(s).', numel(report.unexpectedFailures));
end

end

function lines = section(title, names, results)
lines = {};
if isempty(names)
    return
end
lines{end+1} = '';
lines{end+1} = ['## ' title];
lines{end+1} = '';
for k = 1:numel(names)
    r = results(strcmp({results.Name}, names{k}));
    lines{end+1} = sprintf('- `%s`%s', names{k}, firstDiagnostic(r)); %#ok<AGROW>
end
end

function s = firstDiagnostic(r)
% one line from the failure: the test's own message or the error message
s = '';
if ~(r.Failed || r.Incomplete)
    return
end
try
    rec = r.Details.DiagnosticRecord;
    if isempty(rec)
        return
    end
    rec = rec(1);
    if isprop(rec, 'Exception') && ~isempty(rec.Exception)
        msg = rec.Exception.message;
    elseif ~isempty(rec.TestDiagnosticResults) && strlength(rec.TestDiagnosticResults(1).DiagnosticText) > 0
        msg = rec.TestDiagnosticResults(1).DiagnosticText;
    else
        msg = rec.FrameworkDiagnosticResults(1).DiagnosticText;
    end
    msg = regexprep(char(msg), '\s+', ' ');
    if numel(msg) > 160
        msg = [msg(1:157) '...'];
    end
    s = [' - ' strtrim(msg)];
catch
end
end
