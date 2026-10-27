function [av, info] = mixed_block_sizes()
% Program mixing a 1x1 Gram block (scalar SOS variable s with basis [1]) with a 2x2 matrix SOS
% constraint; regression test for the sparse 4-D indexing in NSD_GRAM (error MATLAB:ndFullOnly).
%   max a s.t. [x^2 + s, a x; a x, 2 + x^2] SOS, s SOS of degree 0, s = 1   ->  a = 1 + sqrt(2)
% The equality keeps the 1x1 Gram block in the interior of its cone.
pvar x a
prog = sosprogram(x, a);
[prog, s] = sossosvar(prog, monomials(x, 0));
prog = soseq(prog, s - 1);
M = poly2dpvar([x^2, a*x; a*x, 2 + x^2], prog.decvartable) + s * [1 0; 0 0];
prog = sosmatrixineq(prog, M, 'Mineq');
prog = sossetobj(prog, -a);
[prog, info] = sossolve_nosdp(prog);
av = double(sosgetsol(prog, a));
end
