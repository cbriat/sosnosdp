function [ok, Phi, info] = nsd_gram(K, E, H, Svec)
%NSD_GRAM  Gram matrices of Theorem 5.2: w = H \ S, Phi_j = Lambda_j(x)^-1 Lambda_j(w) Lambda_j(x)^-1.
%   E = nsd_cone_eval(K, x), H its Hessian. Since sum_j A_j(Phi_j) = H w, the same formula applied
%   to the residual S - sum_j A_j(Phi_j) gives a correction, so the identity is refined iteratively
%   (up to three passes) to rounding level. ok = all Phi_j positive definite and the node identity
%   holds to 1e-8 (relative); callers apply stricter tests (NSD_CHECK_VALUES).
r = K.r; N = K.N;
nb = numel(K.Ps);
Li = cell(1, nb);
for j = 1:nb
  R = E.R{j};
  Li{j} = R \ (R' \ eye(size(R)));
end
Phi = gram_of(K, Li, H \ Svec, r, N);
Srec = apply_A(K, Phi, r, N);
res = norm(Srec - Svec, inf) / max(1e-300, norm(Svec, inf));
for pass = 1:3
  if res < 1e-14, break; end
  D = gram_of(K, Li, H \ (Svec - Srec), r, N);
  Phi_new = cellfun(@(A, B) A + B, Phi, D, 'UniformOutput', false);
  Srec_new = apply_A(K, Phi_new, r, N);
  res_new = norm(Srec_new - Svec, inf) / max(1e-300, norm(Svec, inf));
  if res_new >= res, break; end
  Phi = Phi_new; Srec = Srec_new; res = res_new;
end
mins = zeros(1, nb);
for j = 1:nb
  ev = eig(Phi{j}); mins(j) = min(ev) / max(abs(ev));
end
ok = all(mins > 1e-12) && res < 1e-8;
info = struct('min_eig_rel', min(mins), 'identity_res', res);
end

function Phi = gram_of(K, Li, w, r, N)
Wv = K.Sbig * w; W = reshape(Wv, r, r, N);
[I, J] = ndgrid(1:r, 1:r);
rows = bsxfun(@plus, I(:), (0:N-1)*r); cols = bsxfun(@plus, J(:), (0:N-1)*r);
BD = sparse(rows(:), cols(:), W(:), N*r, N*r);
Phi = cell(1, numel(K.Ps));
for j = 1:numel(K.Ps)
  Pk = K.Pk{j};
  Ph = Li{j} * (Pk' * (BD * Pk)) * Li{j};
  Phi{j} = full((Ph + Ph') / 2);
end
end

function Srec = apply_A(K, Phi, r, N)
rec = zeros(r*r, N);
for j = 1:numel(K.Ps)
  T = full(K.Pk{j} * Phi{j} * K.Pk{j}'); T4 = reshape(T, r, N, r, N);
  for i = 1:N, rec(:, i) = rec(:, i) + reshape(T4(:, i, :, i), r*r, 1); end
end
Srec = reshape(K.S' * rec, [], 1);
end
