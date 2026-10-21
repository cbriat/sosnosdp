function [ok, Phi, info] = nsd_gram(K, E, H, Svec)
%NSD_GRAM  Gram matrices of Theorem 5.2: w = H \ S, Phi_j = Lambda_j(x)^-1 Lambda_j(w) Lambda_j(x)^-1.
%   E = nsd_cone_eval(K, x), H its Hessian. ok = all Phi_j positive definite and the node identity
%   sum_j A_j(Phi_j) = S holds to 1e-8 (relative).
r = K.r; N = K.N;
w = H \ Svec;
Wv = K.Sbig * w; W = reshape(Wv, r, r, N);
[I, J] = ndgrid(1:r, 1:r);
rows = bsxfun(@plus, I(:), (0:N-1)*r); cols = bsxfun(@plus, J(:), (0:N-1)*r);
BD = sparse(rows(:), cols(:), W(:), N*r, N*r);
Phi = cell(1, numel(K.Ps)); mins = zeros(1, numel(K.Ps));
rec = zeros(r*r, N);
for j = 1:numel(K.Ps)
  Pk = K.Pk{j}; R = E.R{j};
  Lw = Pk' * (BD * Pk);
  Li = R \ (R' \ eye(size(R)));
  Ph = Li * Lw * Li; Ph = (Ph + Ph') / 2; Phi{j} = Ph;
  ev = eig(Ph); mins(j) = min(ev) / max(abs(ev));
  T = Pk * Ph * Pk'; T4 = reshape(T, r, N, r, N);
  for i = 1:N, rec(:, i) = rec(:, i) + reshape(T4(:, i, :, i), r*r, 1); end
end
Srec = reshape(K.S' * rec, [], 1);
res = norm(Srec - Svec, inf) / max(1e-300, norm(Svec, inf));
ok = all(mins > 1e-12) && res < 1e-8;
info = struct('min_eig_rel', min(mins), 'identity_res', res);
end
