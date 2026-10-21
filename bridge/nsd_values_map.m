function [cv, Lv] = nsd_values_map(Zact, dimp, c0, C, X, vars, S)
%NSD_VALUES_MAP  Node values of an affine matrix polynomial, stored node by node in svec form.
%   Coefficients are entry-major: index (ent-1)*nZ + k for monomial k and matrix entry ent
%   (column-major over dimp-by-dimp), with value c0 + C*y. Returns cv + Lv*y.
nZ = size(Zact, 1); U = size(X, 1); t = size(S, 2); ny = size(C, 2);
Vz = nsd_moneval(Zact(:, vars), X);                      % U x nZ
V0 = Vz * reshape(full(c0), nZ, dimp^2);                 % U x dimp^2
cv = reshape((V0 * S)', [], 1);
if ny == 0, Lv = zeros(U*t, 0); return; end
Vc = Vz * reshape(full(C), nZ, dimp^2 * ny);             % U x (dimp^2 ny)
Vc = permute(reshape(Vc, U, dimp^2, ny), [2 1 3]);       % dimp^2 x U x ny
Lv = reshape(S' * reshape(Vc, dimp^2, U * ny), t * U, ny);
end
