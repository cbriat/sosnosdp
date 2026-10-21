function sos = quadratic_matrix_ineq(sos, fM)
% What sosmatrixineq(sos, fM) does by default (SOSTOOLS 4, 'quadraticMineq'), with the products
% written as (fM*y).'*y so that they also evaluate in Octave.
n = size(fM, 1);
if isfield(sos, 'varmat') && isfield(sos.varmat, 'vartable') && numel(sos.varmat.vartable) >= n
  names = sos.varmat.vartable(1:n);
else
  names = arrayfun(@(i) sprintf('Mvar_%d', i), (1:n)', 'UniformOutput', false);
  sos.varmat.vartable = names; sos.varmat.symvartable = [];
end
y = polynomial(zeros(n, 1));
for i = 1:n, y(i) = pvar(names{i}); end
sos = sosineq(sos, (fM*y).' * y, 'sparsemultipartite', {sos.vartable, names});
end
