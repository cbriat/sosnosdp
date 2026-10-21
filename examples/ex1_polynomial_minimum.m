%% Example 1: polynomial minimization without SDP
%
% (a) Global minimum of a univariate quartic:   max gam  s.t.  f(x) - gam is SOS.
% (b) Minimum of x^4 - x on the interval [0, 2], imposed directly in the
%     weighted SOS cone (no multiplier, no Gram matrix).

% ---------------------------------------------------------------- (a) global minimum
pvar x gam

prog = sosprogram(x, gam);

f = x^4 - 4*x^3 + 3*x^2 + 2*x + 1;

% f - gam must be a sum of squares
prog = sosineq(prog, f - gam);

% SOSTOOLS minimizes the objective: maximize gam
prog = sossetobj(prog, -gam);

[prog, info] = sossolve_nosdp(prog);

gam_global = double(sosgetsol(prog, gam));

fprintf('global minimum   : %.8f  (%s, %d iterations)\n', gam_global, info.status, info.iters);

% ---------------------------------------------------------------- (b) minimum on [0, 2]
prog = sosprogram(x, gam);

% x^4 - x - gam >= 0 for all x in [0, 2]
prog = sosineq_interval(prog, x^4 - x - gam, x, [0 2]);

prog = sossetobj(prog, -gam);

[prog, info] = sossolve_nosdp(prog);

gam_interval = double(sosgetsol(prog, gam));

% exact value: the minimizer is x = (1/4)^(1/3)
xstar = (1/4)^(1/3);
exact = xstar^4 - xstar;

fprintf('minimum on [0,2] : %.8f  (exact %.8f)\n', gam_interval, exact);
