function [prog, s] = sossosvar_nosdp(prog, Z)
%SOSSOSVAR_NOSDP  SOS decision polynomial without a Gram matrix.
%   [prog, s] = sossosvar_nosdp(prog, Z) declares s = sum of squares of polynomials spanned by the
%   monomial vector Z (as SOSSOSVAR does), but represents s by its coefficients, which are free
%   decision variables, together with the constraint "s is SOS in the basis Z". SOSSOLVE_NOSDP
%   imposes that constraint directly in the interpolant SOS cone, so no Gram matrix appears.
%   The program also remains solvable by SOSTOOLS' own SOSSOLVE.
%   Example:  [prog, s] = sossosvar_nosdp(prog, monomials(x, 0:2));
deg = Z.degmat; names = Z.varname;
if size(deg, 1) ~= numel(Z)
  % Z given as a vector of monomials: one row of degmat per entry
  deg = zeros(numel(Z), numel(names));
  for i = 1:numel(Z)
    zi = Z(i); [~, ix] = ismember(zi.varname, names); deg(i, ix) = full(zi.degmat);
  end
end
deg = full(deg);
nz = size(deg, 1);
[I, J] = ndgrid(1:nz, 1:nz);
ZZ = unique(deg(I(:), :) + deg(J(:), :), 'rows');
mons = polynomial(speye(size(ZZ, 1)), ZZ, names, [size(ZZ, 1) 1]);
[prog, s] = sospolyvar(prog, mons);
prog = sosineq(prog, s);
e = prog.expr.num;
% record the Gram basis in the column order of the program variables
[~, cols] = ismember(names, prog.vartable);
gb = zeros(nz, numel(prog.vartable)); gb(:, cols) = deg;
if ~isfield(prog, 'nosdp') || ~isfield(prog.nosdp, 'gram_basis'), prog.nosdp.gram_basis = {}; end
prog.nosdp.gram_basis{e} = gb;
end
