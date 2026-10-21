function [ok, info] = nsd_check_values(cones, offs, S, delta, EV, mode, depth_tol)
%NSD_CHECK_VALUES  Prove that the node values S lie in the interior of the cone K, with relative depth.
%   EV = cone evaluations at a point x of the dual cone, delta > 0, s_hat = -delta g(x).
%   For each cone k (see NSD_CERTIFY_POINT): Dikin test or explicit Gram matrices. The proved
%   depth is scale-free:
%     Dikin: S - (1 - r) s_hat is in K and s_hat - (delta / lambda_max(Lambda(x))) E is in K, with
%            E = A(I) a fixed interior element, so depth = (1 - r) delta ||E|| / (lambda_max ||S||);
%     Gram : depth = smallest relative eigenvalue of the Gram matrices.
%   mode 'auto' (default) and 'gram': the explicit Gram matrices decide. They must be positive
%   definite with lambda_min > 1e-12 lambda_max, and reproduce S with a relative residual below
%   1e-10 and below 1e-2 times the relative margin lambda_min / lambda_max, so that the residual
%   cannot destroy positive definiteness. This check does not rely on the accuracy of the
%   inverse Hessian. The Dikin ratio is computed and reported, but near the boundary it is not
%   reliable in floating point (a polynomial violating positivity far from the nodes was accepted
%   with depth 1e-8), so in 'auto' mode it does not decide.
%   mode 'dikin': Gram-free; accepted when r < 0.99 and the proved depth is at least depth_tol.
%   Faster, but less robust close to the boundary of the cone.
nk = numel(cones);
ratio = zeros(1, nk); depth = zeros(1, nk); method = cell(1, nk); grams = cell(1, nk);
ok = false;
R = S + delta * EV.g;
HiR = nsd_hsolve(EV, offs, R);
for k = 1:nk
  sl = offs(k)+1:offs(k+1);
  ratio(k) = sqrt(max(0, R(sl)' * HiR(sl))) / delta;
  lam_max = 0;
  for j = 1:numel(EV.E{k}.R), lam_max = max(lam_max, norm(EV.E{k}.R{j})^2); end
  dikin_depth = (1 - ratio(k)) * delta * cones{k}.Enorm / (lam_max * max(norm(S(sl), inf), 1e-300));
  if strcmp(mode, 'dikin') && ratio(k) < 0.99 && dikin_depth >= depth_tol
    method{k} = 'dikin'; depth(k) = dikin_depth;
    continue
  end
  if strcmp(mode, 'dikin')
    info = struct('failed_cone', k, 'dikin_ratio', ratio(k), 'depth', dikin_depth); return;
  end
  [okk, Phi, ginf] = nsd_gram(cones{k}, EV.E{k}, EV.H{k}, S(sl));
  if ~okk || ginf.min_eig_rel < 1e-12 || ginf.identity_res > min(1e-10, 1e-2 * ginf.min_eig_rel)
    info = struct('failed_cone', k, 'dikin_ratio', ratio(k), 'depth', max(dikin_depth, ginf.min_eig_rel), 'detail', ginf);
    return;
  end
  method{k} = 'gram'; grams{k} = Phi; depth(k) = ginf.min_eig_rel;
end
ok = true;
info = struct('gram', {grams}, 'method', {method}, 'dikin_ratio', ratio, 'depth', depth, 'delta', delta);
end
