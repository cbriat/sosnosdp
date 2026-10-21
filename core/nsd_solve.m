function out = nsd_solve(cones, L, c, bobj, opts)
%NSD_SOLVE  Solve  max bobj'z  s.t.  c + L z in K  (K = product of matrix weighted SOS cones).
%   With bobj = 0 (feasibility) the method stops at the first certified iterate.
%   With an objective it runs to convergence and certifies the final iterate.
%   out.status: 'feasible' | 'optimal' | 'infeasible' | 'uncertified' | ...
%   out.z, out.gram (cell per cone of Gram matrices), out.farkas_radius for infeasible problems.
if nargin < 4 || isempty(bobj), bobj = zeros(size(L, 2), 1); end
if nargin < 5, opts = struct(); end
feas = ~any(bobj);
% Exact reparametrization L(:,p) = Q R (column-pivoted QR): the problem depends only on the range
% of L, so solving in w with L z = Q w makes the method independent of the user's basis and
% scaling, and removes redundant variables (Assumption 2.2).
nz = size(L, 2);
[Q, R, p] = qr(full(L), 0);
dR = abs(diag(R));
rank_tol = 1e-14; if isfield(opts, 'rank_tol'), rank_tol = opts.rank_tol; opts = rmfield(opts, 'rank_tol'); end
dmax = max([dR; 1]);
rk = sum(dR > rank_tol * dmax);
near_dep = sum(dR > 1e-16 * dmax & dR < 1e-8 * dmax);   % numerically, not structurally, dependent
Q = Q(:, 1:rk); R11 = R(1:rk, 1:rk); pk = p(1:rk);
bw = R11' \ bobj(pk);
if rk < nz && norm(bobj(p(rk+1:end))) > 1e-9 * max(1, norm(bobj))
  warning('nsd_solve: objective depends on directions that do not affect the constraints');
end
cmode = 'auto';
if isfield(opts, 'certificate'), cmode = opts.certificate; opts = rmfield(opts, 'certificate'); end
depth_tol = 1e-8;
if isfield(opts, 'depth_tol'), depth_tol = opts.depth_tol; opts = rmfield(opts, 'depth_tol'); end
cert = nsd_make_certifier(cones, Q, c, cmode, depth_tol);
offs = [0 cumsum(cellfun(@(K) K.dim, cones))];
o = opts; o.certify = cert; o.stop_on_certificate = feas;
res = nsd_ipm(cones, -Q', bw, c, o);
back = @(w) back_map(w, R11, pk, nz);
out = struct('status', res.status, 'iters', res.iters, 'time', res.time, 'z', [], 'gram', {{}}, ...
             'ipm', rmfield(res, 'EV'), 'cert_method', {{}}, 'dikin_ratio', []);
switch res.status
  case 'certified'
    out.status = 'feasible'; out.z = back(res.info.z); out.gram = res.info.gram;
    out.cert_method = res.info.method; out.dikin_ratio = res.info.dikin_ratio;
  case {'optimal', 'max_iter', 'stalled'}
    % the optimum lies on the boundary of the SOS cone; the certified answer is the last
    % strictly feasible iterate (a rigorous bound), reported with the converged objective
    out.obj_ipm = bw' * (res.y / res.tau);
    if ~isempty(res.last_cert)
      out.z = back(res.last_cert.z); out.gram = res.last_cert.gram; out.cert_iter = res.last_cert.iter;
      out.cert_method = res.last_cert.method; out.dikin_ratio = res.last_cert.dikin_ratio;
      if strcmp(res.status, 'optimal'), out.status = 'optimal'; else, out.status = ['feasible_' res.status]; end
    else
      out.status = ['uncertified_' res.status]; out.z = back(res.y / res.tau);
    end
  case 'infeasible'
    xf = res.x;
    out.farkas_radius = abs(c'*xf) / max(1e-300, norm(Q'*xf));   % radius in the orthonormal coordinates w
end
out.accuracy = res.accuracy;
if ~isempty(out.z) && (~isempty(strfind(out.status, 'feasible')) || strcmp(out.status, 'optimal'))
  % The certificate was obtained for node values computed in the orthonormal coordinates.
  % Re-check it for the values c + L z of the coefficients z actually returned, at the same point.
  % If it fails for the newest certified iterate, fall back to earlier (deeper) ones.
  if strcmp(res.status, 'certified'), cands = {res.info}; else, cands = res.cert_hist; end
  okz = false;
  for ic = numel(cands):-1:1
    cinfo = cands{ic};
    if ~isfield(cinfo, 'x') || isempty(cinfo.x), continue; end
    zc = back(cinfo.z);
    Sz = c + L * zc;
    EVc = nsd_eval_all(cones, offs, cinfo.x);
    [okz, zinfo] = nsd_check_values(cones, offs, Sz, cinfo.delta, EVc, cmode, depth_tol);
    if okz
      out.z = zc; out.gram = cinfo.gram;
      if isfield(cinfo, 'iter'), out.cert_iter = cinfo.iter; end
      out.cert_method = zinfo.method; out.dikin_ratio = zinfo.dikin_ratio; out.depth = zinfo.depth;
      if ~isempty(zinfo.gram) && any(~cellfun(@isempty, zinfo.gram)), out.gram = zinfo.gram; end
      out.coef_consistency = norm(Sz - (c + Q * (Q' * (Sz - c))), inf) / max(1, norm(Sz, inf));
      break
    end
  end
  if ~okz
    out.status = 'uncertified_coefficients';
    out.cert_failure = 'the certificate does not survive the conversion to the returned coefficients';
  end
end
if ~isempty(out.z), out.obj = bobj' * out.z; end
out.rank = rk; out.n_vars = nz; out.near_dependent = near_dep;
if isfield(out, 'obj_ipm') && isfield(out, 'obj'), out.obj_gap = abs(out.obj_ipm - out.obj); end
end

function z = back_map(w, R11, pk, nz)
z = zeros(nz, 1); z(pk) = R11 \ w;
end
