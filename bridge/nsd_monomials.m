function E = nsd_monomials(n, dmin, dmax)
%NSD_MONOMIALS  Exponent rows of all monomials in n variables with total degree in [dmin, dmax].
E = zeros(0, n);
for d = max(dmin, 0):dmax
  E = [E; compositions(n, d)]; %#ok<AGROW>
end
end
function C = compositions(n, d)
if n == 1, C = d; return; end
C = zeros(0, n);
for k = d:-1:0
  R = compositions(n - 1, d - k);
  C = [C; [k * ones(size(R, 1), 1), R]]; %#ok<AGROW>
end
end
