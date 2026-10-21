function prog = sosineq_interval(prog, p, var, ab)
%SOSINEQ_INTERVAL  p(var) >= 0 for var in [a, b], solved directly in the weighted SOS cone
%   with weights 1 and (var - a)(b - var), without multipliers or Gram matrices.
%   Only for programs solved with SOSSOLVE_NOSDP.  p must depend on var only (besides
%   decision variables).  ab = [a b] with a < b finite.
prog = add_interval(prog, p, var, ab);
end
function prog = add_interval(prog, p, var, ab)
assert(numel(ab) == 2 && all(isfinite(ab)) && ab(1) < ab(2), 'interval must be [a b], a < b');
prog = sosconstr(prog, 'ineq', p);
e = prog.expr.num;
prog.expr.type{e} = 'nosdp_interval';
vn = var.varname; if iscell(vn), vn = vn{1}; end
if ~isfield(prog, 'nosdp') || ~isfield(prog.nosdp, 'interval'), prog.nosdp.interval = {}; end
prog.nosdp.interval{e} = struct('var', vn, 'ab', ab(:)');
end
