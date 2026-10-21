%% Example 5: mean-square stabilization of an aperiodic stochastic sampled-data system
%
% Example 16 of C. Briat, Automatica 74 (2016) 279-287:
%
%   dx = (Asd x + Bsd u) dt + Esd x dW1 + alpha Bsd u dW2,   u(t) = Kd [x(t_k); u(t_k)]
%
% with sampling intervals t_{k+1} - t_k in [Tmin, Tmax]. The design conditions of the
% paper's Theorem 15 are two matrix inequalities on intervals (sizes 9 and 6), solved here
% directly in the matrix weighted SOS cone (stochastic_sd_design.m). We
%   (a) design controllers for the ranges of the paper's Table 6,
%   (b) find the largest Tmax for Tmin = 0.001 by bisection,
%   (c) check every design exactly through the second moment E[x x'] (second_moment.m).

% ---------------------------------------------------------------- plant and noise
Asd = [0 1; 0 -1];
Bsd = [0; 1];
Esd = [0 0; 0 0.1];
alpha = 0.1;

n = size(Asd, 1);
mu = size(Bsd, 2);
m = n + mu;

% augmented model (eq. (55) of the paper): state (x, u), u held between samples
sys.Abar = [Asd Bsd; zeros(mu, m)];
sys.E1 = [Esd zeros(n, mu); zeros(mu, m)];
sys.E2 = [zeros(n, n) alpha*Bsd; zeros(mu, m)];
sys.J0 = blkdiag(eye(n), zeros(mu));
sys.Bbar = [zeros(n, mu); eye(mu)];

d = 2;                        % polynomials of order 2, as in the paper
epsl = 1e-3;
x0 = [-1; 1; 0];              % plant state (-1, 1), u(0) = 0

% ---------------------------------------------------------------- (a) ranges of Table 6
ranges = [0.001 0.1;
          0.001 0.5;
          0.001 1;
          1     5;
          1     10;
          1     20];

fprintf('(a) designs for the ranges of Table 6 (degree %d)\n', d);

for i = 1:size(ranges, 1)

  Tmin = ranges(i, 1);
  Tmax = ranges(i, 2);

  [ok, Kd, info] = stochastic_sd_design(sys, Tmin, Tmax, d, epsl);

  if ~ok
    fprintf('  [%g, %g]: %s\n', Tmin, Tmax, info.status);
    continue
  end

  % exact check: 400 random sampling intervals in [Tmin, Tmax]
  rand('seed', i);
  thetas = Tmin + (Tmax - Tmin) * rand(400, 1);
  [msq, rho_p] = second_moment(sys, Kd, x0, thetas);

  % periodic sampling over a grid of the range (necessary condition)
  [~, rho_grid] = second_moment(sys, Kd, x0, linspace(Tmin, Tmax, 50));

  fprintf('  [%5g, %3g]: Kd = %-32s %2d its, %.2f s | E||x||^2: %.1e -> %.1e, max rho = %.4f\n', ...
          Tmin, Tmax, mat2str(Kd, 4), info.iters, info.time, norm(x0)^2, msq(end), max(rho_grid));
end

% ---------------------------------------------------------------- (b) bisection on Tmax
% For the plant above the design stays feasible for very long sampling intervals (the plant is
% marginally stable). To obtain a finite bound we make it unstable, Asd = [0 1; 0 0.5], keep the
% same noise, and search the largest Tmax for Tmin = 0.001 by bisection. Feasibility is monotone
% in Tmax, since a design for [Tmin, Tmax] also works for every smaller Tmax.

sys_u = sys;
sys_u.Abar = [0 1 0; 0 0.5 1; 0 0 0];

Tmin = 0.001;

fprintf('\n(b) unstable variant, largest Tmax for Tmin = %g\n', Tmin);

for d = [2 4]

  lo = 0.5;                   % certified
  hi = 4;                     % not certified

  while (hi - lo) > 1e-3

    mid = (lo + hi) / 2;

    if stochastic_sd_design(sys_u, Tmin, mid, d, epsl)
      lo = mid;
    else
      hi = mid;
    end
  end

  % design at the largest certified Tmax, checked exactly over 20000 random intervals.
  % At the edge of feasibility the stability margin is small: E||x||^2 may first grow
  % before it decays, so a long horizon is needed to see the decay.
  [~, Kd] = stochastic_sd_design(sys_u, Tmin, lo, d, epsl);

  rand('seed', 0);
  thetas = Tmin + (lo - Tmin) * rand(20000, 1);
  msq = second_moment(sys_u, Kd, x0, thetas);

  fprintf('  degree %d: largest certified Tmax = %.4f, Kd = %s\n', d, lo, mat2str(Kd, 4));
  fprintf('            E||x||^2 after 400 / 4000 / 20000 samples: %.1e / %.1e / %.1e\n', ...
          msq(400), msq(4000), msq(end));
end
