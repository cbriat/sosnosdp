%% Example 2: matrix SOS without SDP
%
% Largest a such that
%
%     M(x) = [1 + x^2,   a x  ]
%            [  a x  , 2 + x^2]
%
% is positive semidefinite for all real x. The exact answer is 1 + sqrt(2).
% Both forms of SOSTOOLS' sosmatrixineq are supported by sossolve_nosdp.

pvar x a

M = [1 + x^2, a*x;
     a*x,     2 + x^2];

prog = sosprogram(x, a);

% matrix SOS constraint ('Mineq'); the default form sosmatrixineq(prog, M) also works
prog = sosmatrixineq(prog, M, 'Mineq');

% maximize a
prog = sossetobj(prog, -a);

[prog, info] = sossolve_nosdp(prog);

a_max = double(sosgetsol(prog, a));

fprintf('max a = %.8f  (exact %.8f)\n', a_max, 1 + sqrt(2));
fprintf('cone  : %s\n', info.cones{1});
