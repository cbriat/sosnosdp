function [ok, Kd, info] = stochastic_sd_design(sys, Tmin, Tmax, d, epsl)
%STOCHASTIC_SD_DESIGN  Mean-square stabilization of an aperiodic stochastic sampled-data system.
%
%   Theorem 15 of C. Briat, "Stability analysis and stabilization of stochastic linear impulsive,
%   switched and sampled-data systems under dwell-time constraints", Automatica 74 (2016):
%   find S(tau), tau in [0, Tmax], and Ud such that
%
%     [ S'(tau) + He[Abar S(tau)] ,  S(tau) E1' ,  S(tau) E2' ]
%     [          *                ,  -S(tau)    ,     0       ]  <= 0     for tau in [0, Tmax]   (58)
%     [          *                ,     0       ,  -S(tau)    ]
%
%     [ -S(0) + eps I         ,  (J0 S(0) + Bbar Ud)' ]
%     [ J0 S(0) + Bbar Ud     ,  -S(th)               ]  <= 0         for th in [Tmin, Tmax]  (59)
%
%   Then Kd = Ud S(0)^{-1} makes the sampled-data system mean-square stable for every sequence
%   of sampling intervals in [Tmin, Tmax]. Both conditions are matrix inequalities on intervals,
%   imposed directly in the weighted SOS cone (no multipliers, no SDP).
%
%   sys: struct with fields Abar, E1, E2, J0, Bbar (augmented model of eq. (55)).

m = size(sys.Abar, 1);
mu = size(sys.Bbar, 2);

% normalized clock variable s in [-1, 1] for tau in [0, Tmax]
[s, tm] = nsd_timevar('s', [0 Tmax]);

prog = sosprogram(s);

% decision variables
[prog, S] = sospolymatrixvar(prog, monomials(s, 0:d), [m m], 'symmetric');
[prog, Ud] = sospolymatrixvar(prog, monomials(s, 0), [mu m]);

% ---------------------------------------------------------------- flow condition (58)
% He[Abar S] = Abar S + S Abar' (written with S on the left of the products)
HeAS = (S*sys.Abar')' + S*sys.Abar';

F = [tm.ddt(S) + HeAS,   S*sys.E1',    S*sys.E2';
     (S*sys.E1')',       -S,           zeros(m);
     (S*sys.E2')',       zeros(m),     -S];

% (58) requires F <= 0, i.e. -F >= 0 on the whole clock interval
prog = sosmatrixineq_interval(prog, -F, s, tm.interval());

% ---------------------------------------------------------------- jump condition (59)
S0 = tm.at(S, 0);

% J0 S(0) + Bbar Ud, written as a transpose so that S comes first in the products
W = (S0*sys.J0' + Ud'*sys.Bbar')';

Jm = [-S0 + epsl*eye(m),   W';
      W,                   -S];

% (59) requires Jm <= 0 for every sampling interval th in [Tmin, Tmax]
prog = sosmatrixineq_interval(prog, -Jm, s, tm.interval([Tmin Tmax]));

% ---------------------------------------------------------------- solve
[prog, info] = sossolve_nosdp(prog);

ok = any(strcmp(info.status, {'feasible', 'optimal'}));

Kd = [];
if ok
  S0v = double(sosgetsol(prog, S0));
  Udv = double(sosgetsol(prog, Ud));
  Kd = Udv / S0v;
end
end
