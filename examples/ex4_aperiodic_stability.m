%% Example 4: largest sampling interval under aperiodic sampling, without SDP
%
% x' = A x + B u,  u(t) = K x(t_k) on [t_k, t_{k+1}),  t_{k+1} - t_k in [Tmin, Tmax].
%
% For given Tmax the clock-dependent conditions of aperiodic_feasible.m prove
% stability for every sequence of sampling intervals in [Tmin, Tmax]. We search the
% largest such Tmax by bisection. Feasibility is monotone in Tmax: a certificate on
% [Tmin, Tmax] also covers every smaller Tmax.

% ---------------------------------------------------------------- data (classic benchmark)
A = [0 1; 0 -0.1];
B = [0; 0.1];
K = -[3.75 11.5];

Tmin = 0.01;
epsl = 1e-3;
d = 8;                        % degree of X in the clock variable

% ---------------------------------------------------------------- necessary condition
% periodic sampling with interval th must be stable: rho(G expm(Aa th)) < 1
nx = size(A, 1);
nu = size(B, 2);
m = nx + nu;
Aa = [A B; zeros(nu, m)];
G = blkdiag(eye(nx), zeros(nu)) + [zeros(nx, nu); eye(nu)] * [K zeros(nu)];

th = linspace(Tmin, 3, 3000);
rho = zeros(size(th));
for i = 1:numel(th)
  rho(i) = max(abs(eig(G * expm(Aa * th(i)))));
end
T_periodic = th(find(rho >= 1, 1));

fprintf('periodic sampling loses stability at th = %.4f (upper bound for Tmax)\n', T_periodic);

% ---------------------------------------------------------------- bisection on Tmax
lo = 0.5;                     % certified (checked below)
hi = T_periodic;              % cannot be certified
tol = 1e-3;

if ~aperiodic_feasible(A, B, K, Tmin, lo, d, epsl)
  error('lower bound of the bisection is not certified; decrease lo');
end

while hi - lo > tol

  mid = (lo + hi) / 2;

  [ok, info] = aperiodic_feasible(A, B, K, Tmin, mid, d, epsl);

  fprintf('  Tmax = %.5f : %-10s (%d iterations)\n', mid, info.status, info.iters);

  if ok
    lo = mid;
  else
    hi = mid;
  end
end

fprintf('largest certified Tmax = %.4f  (degree %d, Tmin = %g)\n', lo, d, Tmin);
