function [s, tm] = nsd_timevar(name, ab)
%NSD_TIMEVAR  Normalized independent variable for a physical interval [a, b].
%   [s, tm] = nsd_timevar('s', [a b]) returns a pvar s in [-1, 1] with t = a + (b - a)(s + 1)/2,
%   and helpers to write the program in physical units while SOSTOOLS stores it in s:
%     tm.ddt(Y)          d/dt Y = (2/(b - a)) dY/ds
%     tm.at(Y, t0)       Y evaluated at the physical time t0
%     tm.interval(I)     the s-interval corresponding to a physical interval I (default [a b])
%     tm.t(s0)           physical time of s0
%   Declaring the decision polynomials in s keeps the monomial basis well conditioned, which matters
%   from moderate degrees on (see the README).
%   Example:
%     [s, tm] = nsd_timevar('s', [0 T]);
%     prog = sosprogram(s);
%     [prog, X] = sospolymatrixvar(prog, monomials(s, 0:d), [m m], 'symmetric');
%     prog = sosmatrixineq_interval(prog, tm.ddt(X) - X*Aa' - (X*Aa')', s, tm.interval());
%     XT = tm.at(X, T);  X0 = tm.at(X, 0);
assert(numel(ab) == 2 && ab(1) < ab(2), 'interval must be [a b] with a < b');
a = ab(1); b = ab(2);
s = pvar(name);
tos = @(t0) (2 * t0 - a - b) / (b - a);
tm.a = a; tm.b = b; tm.var = s;
tm.ddt = @(Y) diff(Y, s) * (2 / (b - a));
tm.at = @(Y, t0) subs(Y, s, tos(t0));
tm.interval = @(varargin) interval_of(tos, ab, varargin{:});
tm.t = @(s0) a + (b - a) * (s0 + 1) / 2;
end
function I = interval_of(tos, ab, I)
if nargin < 3, I = ab; end
I = [tos(I(1)), tos(I(2))];
end
