function [ok, info] = aperiodic_feasible(A, B, K, Tmin, Tmax, d, epsl)
%APERIODIC_FEASIBLE  Clock-dependent stability test under aperiodic sampling (used by Example 4).
%
%   ok = true if the conditions below are certified for sampling intervals in [Tmin, Tmax]:
%
%     dX/dtau - Aa X - X Aa' >= 0                          for all tau in [0, Tmax]
%     [X(th), X(th) G'; G X(th), X(0) - eps I] >= 0        for all th  in [Tmin, Tmax]
%
%   with G = J0 + Ba [K 0]. The second constraint is a parameter-dependent LMI.

nx = size(A, 1);
nu = size(B, 2);
m = nx + nu;

% augmented sampled-data model and closed-loop jump map
Aa = [A B; zeros(nu, m)];
J0 = blkdiag(eye(nx), zeros(nu));
Ba = [zeros(nx, nu); eye(nu)];
G = J0 + Ba * [K zeros(nu)];

% normalized clock variable for tau in [0, Tmax]
[s, tm] = nsd_timevar('s', [0 Tmax]);

prog = sosprogram(s);

[prog, X] = sospolymatrixvar(prog, monomials(s, 0:d), [m m], 'symmetric');

% flow condition on [0, Tmax]
P = tm.ddt(X) - X*Aa' - (X*Aa')';
prog = sosmatrixineq_interval(prog, P, s, tm.interval());

% jump condition for every sampling interval th in [Tmin, Tmax]
X0 = tm.at(X, 0);
Psi = [X,         X*G';
       (X*G')',   X0 - epsl*eye(m)];
prog = sosmatrixineq_interval(prog, Psi, s, tm.interval([Tmin Tmax]));

[~, info] = sossolve_nosdp(prog);

ok = any(strcmp(info.status, {'feasible', 'optimal'}));
end
