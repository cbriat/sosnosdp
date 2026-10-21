function [verdict, info] = sampled_data_lmi(A, B, T, d, epsl, opts)
if nargin < 6, opts = struct(); end
% Clock-dependent stabilization LMIs of the paper (Section 6), written in SOSTOOLS:
%   -X'(t) + Aa X + X Aa' <= 0 on [0,T]   (matrix interval constraint, native weighted cone)
%   [X(T), -(J0 X(T) + Ba U)'; -(J0 X(T) + Ba U), X(0) - eps I] >= 0   (constant LMI)
nx = size(A, 1); nu = size(B, 2); m = nx + nu;
Aa = [A B; zeros(nu, m)]; J0 = blkdiag(eye(nx), zeros(nu)); Ba = [zeros(nx, nu); eye(nu)];
pvar t
prog = sosprogram(t);
[prog, X] = sospolymatrixvar(prog, monomials(t, 0:d), [m m], 'symmetric');
[prog, U] = sospolymatrixvar(prog, monomials(t, 0), [nu m]);
Xd = diff(X, t);
P = Xd - X*Aa' - (X*Aa')';
prog = sosmatrixineq_interval(prog, P, t, [0 T]);
XT = subs(X, t, T); X0 = subs(X, t, 0);
Off = -(XT*J0' + U'*Ba')';
Psi = [XT, Off'; Off, X0 - epsl*eye(m)];
prog = sosmatrixineq(prog, Psi, 'Mineq');
[prog, info] = sossolve_nosdp(prog, opts);
verdict = info.status;
if any(strcmp(verdict, {'optimal', 'feasible'})), verdict = 'feasible'; end
if ~isempty(strfind(verdict, 'infeasible')), verdict = 'infeasible'; end
end
