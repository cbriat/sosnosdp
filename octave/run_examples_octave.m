% Run the examples under GNU Octave (see README_OCTAVE.md).
if ~exist('SOSTOOLS_DIR', 'var'), SOSTOOLS_DIR = getenv('SOSTOOLS_DIR'); end
if ~exist('EXAMPLES', 'var'), EXAMPLES = 1:5; end
here = fileparts(mfilename('fullpath'));
addpath(here); addpath(genpath(SOSTOOLS_DIR)); addpath(fullfile(here, '..')); warning('off', 'all');
setup_nosdp();
for k = EXAMPLES
  f = dir(fullfile(here, '..', 'examples', sprintf('ex%d_*.m', k)));
  fprintf('\n===== %s\n', f(1).name); run(fullfile(here, '..', 'examples', f(1).name));
end
