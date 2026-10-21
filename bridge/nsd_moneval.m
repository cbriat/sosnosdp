function V = nsd_moneval(Z, X)
%NSD_MONEVAL  Values of the monomials with exponent rows Z (m-by-n) at the points X (N-by-n).
V = ones(size(X, 1), size(Z, 1));
for v = 1:size(Z, 2)
  e = full(Z(:, v))';
  if any(e), V = V .* bsxfun(@power, X(:, v), e); end
end
end
