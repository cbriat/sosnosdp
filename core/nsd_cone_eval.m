function E = nsd_cone_eval(K, w, need_hess)
%NSD_CONE_EVAL  Membership, barrier gradient and Hessian of the dual cone at w.
%   E.ok     true if every Lambda_j(w) is positive definite
%   E.g      gradient of F(w) = -sum_j log det Lambda_j(w)
%   E.H      Hessian (dense, dim-by-dim), if need_hess
%   E.R{j}   Cholesky factors, Lambda_j = R_j' R_j
if nargin < 3, need_hess = true; end
r = K.r; N = K.N; t = K.t;
E.ok = false; E.g = []; E.H = []; E.R = {}; E.Lam = {};
Wv = K.Sbig * w;                            % stacked vec(W_i), length N r^2
W = reshape(Wv, r, r, N);
% block-diagonal matrix of the W_i
[I, J] = ndgrid(1:r, 1:r);
rows = bsxfun(@plus, I(:), (0:N-1)*r); cols = bsxfun(@plus, J(:), (0:N-1)*r);
BD = sparse(rows(:), cols(:), W(:), N*r, N*r);
g = zeros(K.dim, 1);
if need_hess, H = zeros(K.dim); end
for j = 1:numel(K.Ps)
  Pk = K.Pk{j};
  Lam = Pk' * (BD * Pk); Lam = (Lam + Lam') / 2;
  [R, p] = chol(Lam);
  if p > 0, return; end
  E.R{j} = R; E.Lam{j} = Lam;
  V = R' \ Pk';                             % (L r) x (N r)
  M = V' * V;                               % (N r) x (N r), blocks M_ik
  % gradient: -svec of the diagonal blocks M_ii
  M4 = reshape(M, r, N, r, N);              % M4(p,i,q,k) = M((i-1)r+p, (k-1)r+q)
  Dg = zeros(r*r, N);
  for i = 1:N
    Dg(:, i) = reshape(M4(:, i, :, i), r*r, 1);
  end
  g = g - reshape(K.S' * Dg, [], 1);
  if need_hess
    % block (i,k) of the Hessian in vec coordinates is kron(M_ik, M_ik); then project to svec
    A = bsxfun(@times, reshape(M4, [r 1 N r 1 N]), reshape(M4, [1 r N 1 r N]));
    A = reshape(A, r*r*N, r*r*N);
    H = H + full(K.Sbig' * A * K.Sbig);
  end
end
E.ok = true; E.g = g;
if need_hess, E.H = (H + H') / 2; end
end
