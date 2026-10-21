% scaling_sedumi.m -- clock-dependent sampled-data LMI (oscillator, T = 2) at increasing degree:
%   (1) SeDuMi through sossolve, interval handled by a matrix SOS multiplier (standard SOSTOOLS form)
%   (2) sossolve_nosdp on the same program (multiplier Gram matrix as a PSD cone)
%   (3) sossolve_nosdp with the native interval cone (sosmatrixineq_interval)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'core')); addpath(fullfile(here, '..', 'bridge'));
A = [0.3 3; -3 0.3]; B = [0; 1]; T = 2; epsl = 1e-3;
nx = 2; nu = 1; m = 3;
Aa = [A B; zeros(nu, m)]; J0 = blkdiag(eye(nx), zeros(nu)); Ba = [zeros(nx, nu); eye(nu)];
opt.solver = 'sedumi'; opt.params.fid = 0;
degs = [8 16 24 32];
if exist('scaling_degrees', 'var'), degs = scaling_degrees; end
fprintf('%4s | %-22s | %-22s | %-22s\n', 'd', 'SeDuMi (multiplier)', 'no SDP (multiplier)', 'no SDP (native)');
for d = degs
  row = cell(1, 3);
  for mode = 1:3
    pvar t
    prog = sosprogram(t);
    [prog, X] = sospolymatrixvar(prog, monomials(t, 0:d), [m m], 'symmetric');
    [prog, U] = sospolymatrixvar(prog, monomials(t, 0), [nu m]);
    P = diff(X, t) - X*Aa' - (X*Aa')';
    if mode == 3
      prog = sosmatrixineq_interval(prog, P, t, [0 T]);
    else
      [prog, Sm] = sossosmatrixvar(prog, monomials(t, 0:d/2-1), [m m]);
      prog = sosmatrixineq(prog, P - Sm*(t*(T - t)), 'Mineq');
    end
    XT = subs(X, t, T); X0 = subs(X, t, 0);
    Off = -(XT*J0' + U'*Ba')';
    prog = sosmatrixineq(prog, [XT, Off'; Off, X0 - epsl*eye(m)], 'Mineq');
    t0 = tic;
    if mode == 1
      prog = sossolve(prog, opt); el = toc(t0);
      st = prog.solinfo.info; ok = st.numerr == 0 && st.pinf == 0 && st.dinf == 0;
      row{mode} = sprintf('%-9s %8.2f s', ternary(ok, 'feasible', sprintf('numerr=%d pinf=%d dinf=%d', st.numerr, st.pinf, st.dinf)), el);
    else
      [prog, info] = sossolve_nosdp(prog); el = toc(t0);
      row{mode} = sprintf('%-9s %8.2f s', info.status, el);
    end
  end
  fprintf('%4d | %-22s | %-22s | %-22s\n', d, row{:});
end
