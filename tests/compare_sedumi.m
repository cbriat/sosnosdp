% compare_sedumi.m -- the same SOSTOOLS programs solved by sossolve (SeDuMi) and sossolve_nosdp.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'core')); addpath(fullfile(here, '..', 'bridge'));
opt.solver = 'sedumi'; opt.params.fid = 0;
fprintf('%-44s %14s %14s %9s %9s\n', 'problem', 'SeDuMi', 'no SDP', 't SeDuMi', 't no SDP');
pvar x gam a
probs = {'univariate minimum', 'matrix SOS (Mineq), max a', 'half-line [1,inf)', 'interval [0,2] via sossosvar multiplier'};
for k = 1:numel(probs)
  vals = zeros(1, 2); ts = zeros(1, 2);
  for mode = 1:2
    switch k
      case 1
        prog = sosprogram(x, gam); prog = sosineq(prog, x^4 - 4*x^3 + 3*x^2 + 2*x + 1 - gam); prog = sossetobj(prog, -gam); V = gam;
      case 2
        prog = sosprogram(x, a); prog = sosmatrixineq(prog, [1 + x^2, a*x; a*x, 2 + x^2], 'Mineq'); prog = sossetobj(prog, -a); V = a;
      case 3
        prog = sosprogram(x, gam); prog = sosineq(prog, x^3 - 3*x - gam, [1 inf]); prog = sossetobj(prog, -gam); V = gam;
      case 4
        prog = sosprogram(x, gam); [prog, s] = sossosvar(prog, [1; x]);
        prog = sosineq(prog, poly2dpvar(x^4 - x - gam, prog.decvartable) - s * (x*(2 - x))); prog = sossetobj(prog, -gam); V = gam;
    end
    t0 = tic;
    if mode == 1, prog = sossolve(prog, opt); else, prog = sossolve_nosdp(prog); end
    ts(mode) = toc(t0); vals(mode) = double(sosgetsol(prog, V));
  end
  fprintf('%-44s %14.8f %14.8f %8.2fs %8.2fs\n', probs{k}, vals(1), vals(2), ts(1), ts(2));
end
