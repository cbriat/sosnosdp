function S = nsd_svec_basis(r)
%NSD_SVEC_BASIS  Orthonormal basis of symmetric r-by-r matrices, as columns of vec(E_a).
%   svec(X) = S' * X(:) and X(:) = S * svec(X). The order is the upper triangle column by
%   column, (1,1),(1,2),(2,2),(1,3),..., and off-diagonal entries carry the factor sqrt(2).
t = r*(r+1)/2;
S = zeros(r*r, t);
a = 0;
for q = 1:r
  for p = 1:q
    a = a + 1;
    if p == q
      S((q-1)*r + p, a) = 1;
    else
      S((q-1)*r + p, a) = 1/sqrt(2);
      S((p-1)*r + q, a) = 1/sqrt(2);
    end
  end
end
end
