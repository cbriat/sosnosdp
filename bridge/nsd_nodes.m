function X = nsd_nodes(ZZ, lo, hi)
%NSD_NODES  Unisolvent nodes for span{monomials ZZ} (U = rows of ZZ) in the box [lo, hi]^n:
%   approximate Fekete points (Sommariva-Vianello), selected by column-pivoted QR from a large
%   candidate set (Chebyshev-Lobatto grid for one variable, Halton points otherwise).
[U, n] = size(ZZ);
if nargin < 2, lo = -ones(1, n); hi = ones(1, n); end
if n == 1
  m = max(40 * U, 400);
  cand = cos(pi * (0:m-1)' / (m - 1));
else
  m = max(30 * U, 1000);
  cand = 2 * nsd_halton(m, n) - 1;
end
cand = bsxfun(@plus, (lo + hi) / 2, bsxfun(@times, (hi - lo) / 2, cand));
V = nsd_moneval(ZZ, cand);
sc = max(abs(V), [], 1); sc(sc == 0) = 1;
[~, ~, p] = qr(bsxfun(@rdivide, V, sc)', 0);
X = cand(p(1:U), :);
Vs = nsd_moneval(ZZ, X);
if rcond(Vs) < 1e-14
  warning('nsd_nodes: interpolation matrix is ill-conditioned (rcond %.1e)', rcond(Vs));
end
end
