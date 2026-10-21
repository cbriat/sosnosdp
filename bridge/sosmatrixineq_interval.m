function prog = sosmatrixineq_interval(prog, M, var, ab)
%SOSMATRIXINEQ_INTERVAL  M(var) positive semidefinite for var in [a, b], solved directly in the
%   matrix weighted SOS cone (weights 1 and (var - a)(b - var)), without multipliers or Gram
%   matrices. Only for programs solved with SOSSOLVE_NOSDP.
assert(size(M, 1) == size(M, 2), 'M must be square');
prog = sosineq_interval(prog, M, var, ab);
end
