% Run the test suite under GNU Octave. Set SOSTOOLS_DIR to a SOSTOOLS 4 copy patched as described
% in README_OCTAVE.md. MATLAB users run tests/run_tests.m directly.
if ~exist('SOSTOOLS_DIR', 'var'), SOSTOOLS_DIR = getenv('SOSTOOLS_DIR'); end
here = fileparts(mfilename('fullpath'));
addpath(here); addpath(genpath(SOSTOOLS_DIR)); warning('off', 'all');
run(fullfile(here, '..', 'tests', 'run_tests.m'));
