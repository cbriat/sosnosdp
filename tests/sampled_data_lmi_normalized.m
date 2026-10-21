function [verdict, info, prog, X, U] = sampled_data_lmi_normalized(A, B, T, d, epsl, opts)
% Clock-dependent stabilization LMIs with normalized time s in [-1, 1], t = T (s + 1) / 2:
%   -(2/T) X'(s) + Aa X + X Aa' <= 0 on [-1, 1],  boundary LMI at s = 1 and s = -1.
if nargin < 6, opts = struct(); end
nx = size(A, 1); nu = size(B, 2); m = nx + nu;
Aa = [A B; zeros(nu, m)]; J0 = blkdiag(eye(nx), zeros(nu)); Ba = [zeros(nx, nu); eye(nu)];
pvar s
prog = sosprogram(s);
[prog, X] = sospolymatrixvar(prog, monomials(s, 0:d), [m m], 'symmetric');
[prog, U] = sospolymatrixvar(prog, monomials(s, 0), [nu m]);
P = diff(X, s) * (2 / T) - X*Aa' - (X*Aa')';
prog = sosmatrixineq_interval(prog, P, s, [-1 1]);
XT = subs(X, s, 1); X0 = subs(X, s, -1);
Off = -(XT*J0' + U'*Ba')';
prog = sosmatrixineq(prog, [XT, Off'; Off, X0 - epsl*eye(m)], 'Mineq');
[prog, info] = sossolve_nosdp(prog, opts);
verdict = info.status;
if any(strcmp(verdict, {'optimal', 'feasible'})), verdict = 'feasible'; end
if ~isempty(strfind(verdict, 'infeasible')), verdict = 'infeasible'; end
end
