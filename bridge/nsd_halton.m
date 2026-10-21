function X = nsd_halton(m, n)
%NSD_HALTON  First m points of the n-dimensional Halton sequence in (0,1)^n (deterministic).
p = primes(200); p = p(1:n);
X = zeros(m, n);
for j = 1:n
  for i = 1:m
    f = 1; r = 0; k = i;
    while k > 0
      f = f / p(j); r = r + f * mod(k, p(j)); k = floor(k / p(j));
    end
    X(i, j) = r;
  end
end
end
