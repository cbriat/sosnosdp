% run_tests.m -- tests of sossolve_nosdp against problems with known answers.
% Products are written with dpvar objects on the left, which also works in Octave.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..')); setup_nosdp();
npass = 0; nfail = 0;
report = @(name, ok, msg) fprintf('%-58s %s  %s\n', name, ternary(ok, 'PASS', 'FAIL'), msg);

% ---- T1: univariate minimum (exact: nonnegative univariate polynomials are SOS)
pvar x gam
prog = sosprogram(x, gam);
f = x^4 - 4*x^3 + 3*x^2 + 2*x + 1;
prog = sosineq(prog, f - gam); prog = sossetobj(prog, -gam);
[prog, info] = sossolve_nosdp(prog);
g = double(sosgetsol(prog, gam));
r = roots([4 -12 6 2]); r = r(abs(imag(r)) < 1e-9); fmin = min(polyval([1 -4 3 2 1], real(r)));
% certified bound: the last iterate with a proved depth inside the cone, slightly below the optimum
ok = strcmp(info.status, 'optimal') && g <= fmin + 1e-9 && fmin - g < 1e-5;
report('T1 univariate minimum (sosineq, sosgetsol)', ok, sprintf('gamma = %.8f, true min = %.8f, %d its', g, fmin, info.iters));
npass = npass + ok; nfail = nfail + ~ok;

% ---- T2: unbounded polynomial -> infeasible
pvar x1 x2
prog = sosprogram([x1; x2], gam);
prog = sosineq(prog, x1^4 + x2^4 - 3*x1^2*x2^2 + x1^2 + 2 - gam); prog = sossetobj(prog, -gam);
[prog, info] = sossolve_nosdp(prog);
ok = ~isempty(strfind(info.status, 'infeasible')) || ~isempty(strfind(info.status, 'uncertified'));
report('T2 unbounded quartic (must not be certified)', ok, sprintf('status = %s', info.status));
npass = npass + ok; nfail = nfail + ~ok;

% ---- T3: matrix SOS with 'Mineq': max a s.t. [1+x^2, a x; a x, 2+x^2] SOS  ->  1 + sqrt(2)
pvar x a
prog = sosprogram(x, a);
M = [1 + x^2, a*x; a*x, 2 + x^2];
prog = sosmatrixineq(prog, M, 'Mineq'); prog = sossetobj(prog, -a);
[prog, info] = sossolve_nosdp(prog);
av = double(sosgetsol(prog, a));
ok = strcmp(info.status, 'optimal') && abs(av - (1 + sqrt(2))) < 1e-6;
report('T3 matrix SOS, Mineq (max a = 1+sqrt(2))', ok, sprintf('a = %.8f, %d its', av, info.iters));
npass = npass + ok; nfail = nfail + ~ok;

% ---- T4: same problem in SOSTOOLS' default form y'M(x)y (sparsemultipartite, auxiliary vars)
prog = sosprogram(x, a);
if exist('OCTAVE_VERSION', 'builtin')
  prog = quadratic_matrix_ineq(prog, M);   % Octave: equivalent reconstruction (see octave/README_OCTAVE.md)
else
  prog = sosmatrixineq(prog, M);           % MATLAB: SOSTOOLS' default form
end
prog = sossetobj(prog, -a);
[prog, info] = sossolve_nosdp(prog);
av = double(sosgetsol(prog, a));
ok = strcmp(info.status, 'optimal') && abs(av - (1 + sqrt(2))) < 1e-6;
report('T4 matrix SOS, quadraticMineq form (y''My)', ok, sprintf('a = %.8f, cone: %s', av, info.cones{1}));
npass = npass + ok; nfail = nfail + ~ok;

% ---- T5: half-line: max gam s.t. x^3 - 3x - gam >= 0 on [1, inf)  ->  -2
pvar x gam
prog = sosprogram(x, gam);
prog = sosineq(prog, x^3 - 3*x - gam, [1 inf]); prog = sossetobj(prog, -gam);
[prog, info] = sossolve_nosdp(prog);
g = double(sosgetsol(prog, gam));
ok = strcmp(info.status, 'optimal') && abs(g + 2) < 1e-6;
report('T5 half-line constraint (SOSTOOLS posint)', ok, sprintf('gamma = %.8f, %s', g, info.status));
npass = npass + ok; nfail = nfail + ~ok;

% ---- T6a: finite interval, native weighted cone: min x^4 - x on [0,2]
xs = (1/4)^(1/3); fmin = xs^4 - xs;
prog = sosprogram(x, gam);
prog = sosineq_interval(prog, x^4 - x - gam, x, [0 2]); prog = sossetobj(prog, -gam);
[prog, info] = sossolve_nosdp(prog);
g = double(sosgetsol(prog, gam));
ok = strcmp(info.status, 'optimal') && abs(g - fmin) < 1e-6;
report('T6a interval [0,2], native (sosineq_interval)', ok, sprintf('gamma = %.8f, true = %.8f', g, fmin));
npass = npass + ok; nfail = nfail + ~ok;

% ---- T6b: same via a Putinar multiplier s(x) SOS:  f - gam - s x (2 - x) SOS
prog = sosprogram(x, gam);
[prog, s] = sossosvar(prog, [1; x]);
prog = sosineq(prog, poly2dpvar(x^4 - x - gam, prog.decvartable) - s * (x*(2 - x)));
prog = sossetobj(prog, -gam);
[prog, info] = sossolve_nosdp(prog);
g = double(sosgetsol(prog, gam));
ok = strcmp(info.status, 'optimal') && abs(g - fmin) < 1e-6;
report('T6b interval via sossosvar multiplier', ok, sprintf('gamma = %.8f, cones: %d, %s', g, numel(info.cones), info.status));
npass = npass + ok; nfail = nfail + ~ok;

% ---- T6c: equality constraint: same as T6b but the multiplier is defined through soseq
prog = sosprogram(x, gam);
[prog, s] = sossosvar(prog, [1; x]);
[prog, q] = sospolyvar(prog, monomials(x, 0:4));
prog = soseq(prog, q - (poly2dpvar(x^4 - x - gam, prog.decvartable) - s * (x*(2 - x))));
prog = sosineq(prog, q); prog = sossetobj(prog, -gam);
[prog, info] = sossolve_nosdp(prog);
g = double(sosgetsol(prog, gam));
ok = strcmp(info.status, 'optimal') && abs(g - fmin) < 1e-6;
report('T6c with soseq and sospolyvar', ok, sprintf('gamma = %.8f, %s', g, info.status));
npass = npass + ok; nfail = nfail + ~ok;

% ---- T7: clock-dependent sampled-data LMI (oscillator, T = 2): feasible at d = 8, not at d = 6
A = [0.3 3; -3 0.3]; B = [0; 1]; T = 2; epsl = 1e-3;
for d = [6 8]
  [verdict, inf7] = sampled_data_lmi(A, B, T, d, epsl);
  expect = ternary(d == 8, 'feasible', 'infeasible');
  ok = strcmp(verdict, expect);
  report(sprintf('T7 clock-dependent LMI, oscillator T=2, d=%d', d), ok, sprintf('%s (expected %s), %d its, %.2f s', inf7.status, expect, inf7.iters, inf7.time));
  npass = npass + ok; nfail = nfail + ~ok;
end
% ---- T8: SOS multiplier without Gram matrix (sossosvar_nosdp), interval via Putinar form
pvar x gam
prog = sosprogram(x, gam);
[prog, s8] = sossosvar_nosdp(prog, [1; x]);
prog = sosineq(prog, poly2dpvar(x^4 - x - gam, prog.decvartable) - s8 * (x*(2 - x)));
prog = sossetobj(prog, -gam);
[prog, info] = sossolve_nosdp(prog);
g = double(sosgetsol(prog, gam));
psd_cones = sum(~cellfun(@isempty, strfind(info.cones, 'Gram')));
ok = strcmp(info.status, 'optimal') && abs(g - fmin) < 1e-6 && psd_cones == 0;
report('T8 multiplier via sossosvar_nosdp (no Gram matrix)', ok, sprintf('gamma = %.8f, Gram cones: %d', g, psd_cones));
npass = npass + ok; nfail = nfail + ~ok;

% ---- T9: clock-dependent LMI at high degree with normalized time s in [-1, 1]
for d = [24 40]
  [verdict, inf9] = sampled_data_lmi_normalized(A, B, T, d, epsl);
  ok = strcmp(verdict, 'feasible');
  report(sprintf('T9 clock-dependent LMI, normalized time, d=%d', d), ok, sprintf('%s, %d its, %.2f s, certificate: %s', inf9.status, inf9.iters, inf9.time, strjoin(unique(inf9.cert_method), '+')));
  npass = npass + ok; nfail = nfail + ~ok;
end

% ---- T10: Gram-free certificate only ('dikin') on the d = 8 instance
[verdict, inf10] = sampled_data_lmi(A, B, T, 8, epsl, struct('certificate', 'dikin'));
ok = strcmp(verdict, 'feasible') && all(strcmp(inf10.cert_method, 'dikin'));
report('T10 Gram-free certificate only (Dikin test)', ok, sprintf('%s, %d its, Dikin ratios %s', inf10.status, inf10.iters, mat2str(round(inf10.dikin_ratio*100)/100)));
npass = npass + ok; nfail = nfail + ~ok;

% ---- T11: nsd_timevar at high degree, physical-time formulation; check the recovered controller
nxx = size(A, 1); m = nxx + size(B, 2);
Aa = [A B; zeros(size(B, 2), m)]; J0 = blkdiag(eye(nxx), zeros(size(B, 2))); Ba = [zeros(nxx, size(B, 2)); eye(size(B, 2))];
[s, tm] = nsd_timevar('s', [0 T]);
prog = sosprogram(s);
[prog, X] = sospolymatrixvar(prog, monomials(s, 0:32), [m m], 'symmetric');
[prog, U] = sospolymatrixvar(prog, monomials(s, 0), [size(B, 2) m]);
prog = sosmatrixineq_interval(prog, tm.ddt(X) - X*Aa' - (X*Aa')', s, tm.interval());
XT = tm.at(X, T); X0 = tm.at(X, 0); Off = -(XT*J0' + U'*Ba')';
prog = sosmatrixineq(prog, [XT, Off'; Off, X0 - epsl*eye(m)], 'Mineq');
[prog, info] = sossolve_nosdp(prog);
K = double(sosgetsol(prog, U)) / double(sosgetsol(prog, XT));
rho = max(abs(eig((J0 + Ba*K) * expm(Aa*T))));
ok = strcmp(info.status, 'feasible') && info.rank == info.n_free && rho < 1;
report('T11 nsd_timevar, d=32, recovered controller', ok, sprintf('%s, rank %d of %d, rho = %.4f', info.status, info.rank, info.n_free, rho));
npass = npass + ok; nfail = nfail + ~ok;

% ---- T12: regression test for a false positive found earlier: with eps = 0 the d = 6 problem
% (infeasible with any eps > 0) admits only solutions on the boundary of the cone; a near-boundary
% point must not be accepted as a certificate.
[s, tm] = nsd_timevar('s', [0 T]);
prog = sosprogram(s);
[prog, X] = sospolymatrixvar(prog, monomials(s, 0:6), [m m], 'symmetric');
[prog, U] = sospolymatrixvar(prog, monomials(s, 0), [size(B, 2) m]);
prog = sosmatrixineq_interval(prog, tm.ddt(X) - X*Aa' - (X*Aa')', s, tm.interval());
XT = tm.at(X, T); X0 = tm.at(X, 0); Off = -(XT*J0' + U'*Ba')';
prog = sosmatrixineq(prog, [XT, Off'; Off, X0], 'Mineq');
[prog, info] = sossolve_nosdp(prog);
ok = ~any(strcmp(info.status, {'feasible', 'optimal'}));
report('T12 no certificate on the boundary (eps = 0, d = 6)', ok, sprintf('status = %s', info.status));
npass = npass + ok; nfail = nfail + ~ok;

% ---- T13: global matrix SOS in two variables, max a with exact value 1 + 1/sqrt(2)
pvar x y z a
M = [x^2 + y^2 + 1, a*(x + y); a*(x + y), 2 + x^2 + y^2];
prog = sosprogram([x; y], a);
prog = sosmatrixineq(prog, M, 'Mineq');
prog = sossetobj(prog, -a);
[prog, info] = sossolve_nosdp(prog);
av = double(sosgetsol(prog, a)); ex = 1 + 1/sqrt(2);
ok = strcmp(info.status, 'optimal') && av <= ex + 1e-9 && ex - av < 1e-4;
report('T13 matrix SOS, 2 variables (exact 1+1/sqrt(2))', ok, sprintf('a = %.8f, exact %.8f, %s', av, ex, info.cones{1}));
npass = npass + ok; nfail = nfail + ~ok;

% ---- T14: global scalar SOS in three variables, x^4+y^4+z^4+1+a xyz, max a = 4
prog = sosprogram([x; y; z], a);
prog = sosineq(prog, x^4 + y^4 + z^4 + 1 + a*x*y*z);
prog = sossetobj(prog, -a);
[prog, info] = sossolve_nosdp(prog);
av = double(sosgetsol(prog, a));
ok = strcmp(info.status, 'optimal') && av <= 4 + 1e-9 && 4 - av < 1e-4;
report('T14 scalar SOS, 3 variables (exact 4)', ok, sprintf('a = %.8f, %s', av, info.cones{1}));
npass = npass + ok; nfail = nfail + ~ok;

% ---- T15: a matrix polynomial on the boundary of the SOS cone for every a (constant (2,2)
% entry forces a singular Gram matrix); no interior certificate exists and none may be claimed.
% max a would be 1 (the matrix is not PSD for a > 1).
M = [x^2 + y^2 + 1, a*(x + y); a*(x + y), 2];
prog = sosprogram([x; y], a);
prog = sosmatrixineq(prog, M, 'Mineq');
prog = sossetobj(prog, -a);
[prog, info] = sossolve_nosdp(prog);
certified = any(strcmp(info.status, {'feasible', 'optimal'}));
av = double(sosgetsol(prog, a));
ok = ~certified || av <= 1 + 1e-9;
report('T15 structurally singular Gram matrix: no false certificate', ok, sprintf('status = %s, a = %.7f', info.status, av));
npass = npass + ok; nfail = nfail + ~ok;

% ---- T16: box in two variables through Gram-free multipliers: min xy on [-1,1]^2 = -1
pvar gam
prog = sosprogram([x; y], gam);
Zm = monomials([x; y], 0:1);
[prog, s1] = sossosvar_nosdp(prog, Zm);
[prog, s2] = sossosvar_nosdp(prog, Zm);
prog = sosineq(prog, poly2dpvar(x*y - gam, prog.decvartable) - s1*(1 - x^2) - s2*(1 - y^2));
prog = sossetobj(prog, -gam);
[prog, info] = sossolve_nosdp(prog);
g = double(sosgetsol(prog, gam));
ok = strcmp(info.status, 'optimal') && g <= -1 + 1e-9 && g > -1 - 1e-5;
report('T16 box [-1,1]^2 via sossosvar_nosdp multipliers', ok, sprintf('gamma = %.8f (exact -1), %d cones', g, numel(info.cones)));
npass = npass + ok; nfail = nfail + ~ok;

% ---- T17: matrix inequality on a box in two variables through matrix SOS multipliers
% [1+x^2, a x y; a x y, 1+y^2] >= 0 on [-1,1]^2 iff a <= 2
prog = sosprogram([x; y], a);
[prog, S1] = sossosmatrixvar(prog, Zm, [2 2]);
[prog, S2] = sossosmatrixvar(prog, Zm, [2 2]);
M = poly2dpvar([1 + x^2, a*x*y; a*x*y, 1 + y^2], prog.decvartable);
prog = sosmatrixineq(prog, M - S1*(1 - x^2) - S2*(1 - y^2), 'Mineq');
prog = sossetobj(prog, -a);
[prog, info] = sossolve_nosdp(prog);
av = double(sosgetsol(prog, a));
ok = strcmp(info.status, 'optimal') && av <= 2 + 1e-9 && 2 - av < 1e-5;
report('T17 matrix inequality on [-1,1]^2 via sossosmatrixvar', ok, sprintf('a = %.8f (exact 2), %d cones', av, numel(info.cones)));
npass = npass + ok; nfail = nfail + ~ok;

% ---- T18: regression test, a cone of dimension one (here the 1x1 Gram block of a degree-0
% sossosvar, next to a 2x2 matrix SOS constraint) made NSD_GRAM index a sparse matrix as a 4-D
% array (error MATLAB:ndFullOnly)
[av, info] = mixed_block_sizes();
ok = strcmp(info.status, 'optimal') && abs(av - (1 + sqrt(2))) < 1e-6;
report('T18 mixed Gram block sizes (1x1 and 2x2)', ok, sprintf('a = %.8f (exact 1+sqrt(2)), %d cones', av, numel(info.cones)));
npass = npass + ok; nfail = nfail + ~ok;

fprintf('\n%d passed, %d failed\n', npass, nfail);
if nfail > 0
  error('sosnosdp:tests', '%d test(s) failed', nfail);
end
