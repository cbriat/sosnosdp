function setup_nosdp()
%SETUP_NOSDP  Add SOSNOSDP to the path and check the SOSTOOLS installation.
%   Run once per session after adding SOSTOOLS:  addpath(genpath('<SOSTOOLS folder>')); setup_nosdp
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'core'), fullfile(root, 'bridge'), fullfile(root, 'examples'));
if exist('sosprogram', 'file') ~= 2
  warning('setup_nosdp: SOSTOOLS was not found. Run addpath(genpath(''<SOSTOOLS folder>'')) first.');
  return
end
% SOSTOOLS' sosprogram refuses to build a program unless an SDP solver is on the path.
solvers = {'sedumi', 'sqlp', 'csdp', 'sdpnalplus', 'sdpnal', 'cdcs', 'sdpam', 'mosekopt'};
if ~any(cellfun(@(f) exist(f, 'file') > 0, solvers))
  addpath(fullfile(root, 'placeholder'));
  fprintf('No SDP solver found: added a placeholder so that SOSTOOLS can build programs.\n');
end
fprintf('SOSNOSDP ready. Solve programs with  [prog, info] = sossolve_nosdp(prog).\n');
end
