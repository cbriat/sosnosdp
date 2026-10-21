function res = nsd_ipm(cones, A, b, c, opts)
%NSD_IPM  Homogeneous predictor-corrector interior-point method (Skajaa-Ye, as in Papp-Yildiz)
%   primal  min c'x  s.t. A x = b,  x in K*      (K* carries the barrier -sum log det Lambda_j)
%   dual    max b'y  s.t. s = c - A'y in K       (K = product of matrix weighted SOS cones)
%   opts: max_iter (200), tol (1e-8), eta_pred (0.7), eta_corr (0.3), max_corr (6),
%         certify (function handle @(x,y,tau,kappa,mu,EV) -> [ok, info]), verbose (false)
if nargin < 5, opts = struct(); end
o = struct('max_iter', 200, 'tol', 1e-8, 'eta_pred', 0.70, 'eta_corr', 0.30, 'max_corr', 6, ...
           'certify', [], 'stop_on_certificate', true, 'verbose', false);
f = fieldnames(opts); for i = 1:numel(f), o.(f{i}) = opts.(f{i}); end
t0 = tic;
offs = [0 cumsum(cellfun(@(K) K.dim, cones))];
nu = sum(cellfun(@(K) K.nu, cones));
my = size(A, 1);
x = cell2mat(cellfun(@(K) K.x0, cones(:), 'UniformOutput', false));
EV = nsd_eval_all(cones, offs, x);
s = -EV.g; y = zeros(my, 1); tau = 1; kappa = 1;
muf = @(x, s, tau, kappa) (x'*s + tau*kappa) / (nu + 1);
rp0 = norm(A*x - b*tau) + 1; rd0 = norm(-A'*y + c*tau - s) + 1;
status = 'max_iter'; info = struct(); hist = zeros(0, 6); last_cert = []; last_fail = []; cert_hist = {}; accuracy = '';
for it = 1:o.max_iter
  mu = muf(x, s, tau, kappa);
  rp = A*x - b*tau; rd = -A'*y + c*tau - s; rg = b'*y - c'*x - kappa;
  pr = prox(EV, offs, s, tau, kappa, mu);
  hist(end+1, :) = [it, mu, tau, kappa, norm(rp)/rp0, norm(rd)/rd0]; %#ok<AGROW>
  if o.verbose
    fprintf('%3d mu=%.2e tau=%.2e kap=%.2e rp=%.1e rd=%.1e prox=%.2f\n', it, mu, tau, kappa, norm(rp)/rp0, norm(rd)/rd0, pr);
  end
  if ~isempty(o.certify) && tau > kappa
    [okc, infoc] = o.certify(x, y, tau, kappa, mu, EV);
    if okc
      if o.stop_on_certificate, status = 'certified'; info = infoc; break; end
      last_cert = infoc; last_cert.iter = it;
      cert_hist{end+1} = last_cert;
      if numel(cert_hist) > 10, cert_hist(1) = []; end
    else
      last_fail = infoc; last_fail.iter = it;
    end
  end
  if mu < o.tol && tau < o.tol * max(1, kappa) && c'*x < 0
    status = 'infeasible'; info = struct('cx', c'*x, 'Ax', norm(A*x)); break;
  end
  gapr = mu / max(tau, 1e-300)^2; resr = max(norm(rp)/rp0, norm(rd)/rd0);
  if resr < o.tol && gapr < o.tol && tau > kappa
    status = 'optimal'; accuracy = 'full'; break;
  end
  if resr < 1e3 * o.tol && gapr < o.tol && tau > kappa
    % residuals stagnate near the boundary while the gap keeps shrinking (loss of accuracy in
    % the Hessian); stop with reduced accuracy, as SDP solvers do
    status = 'optimal'; accuracy = 'reduced'; break;
  end
  % predictor
  [dx, dy, dtau, ds, dkap] = direction(EV, offs, A, b, c, mu, tau, -rp, -rd, -rg, -s, -kappa);
  alpha = 1;
  if dtau < 0, alpha = min(alpha, 0.99 * tau / -dtau); end
  if dkap < 0, alpha = min(alpha, 0.99 * kappa / -dkap); end
  accepted = false;
  for ls = 1:30
    xn = x + alpha*dx;
    EVn = nsd_eval_all(cones, offs, xn);
    if EVn.ok
      sn = s + alpha*ds; tn = tau + alpha*dtau; kn = kappa + alpha*dkap;
      mun = muf(xn, sn, tn, kn);
      if mun > 0 && prox(EVn, offs, sn, tn, kn, mun) <= o.eta_pred
        accepted = true; break;
      end
    end
    alpha = 0.7 * alpha;
  end
  if ~accepted, status = 'stalled'; break; end
  x = xn; y = y + alpha*dy; s = sn; tau = tn; kappa = kn; EV = EVn;
  % correctors
  for kc = 1:o.max_corr
    mu = muf(x, s, tau, kappa);
    if prox(EV, offs, s, tau, kappa, mu) <= o.eta_corr, break; end
    [dx, dy, dtau, ds, dkap] = direction(EV, offs, A, b, c, mu, tau, 0*rp, 0*rd, 0, -s - mu*EV.g, -kappa + mu/tau);
    beta = 1; okstep = false;
    for ls = 1:20
      xn = x + beta*dx; tn = tau + beta*dtau; kn = kappa + beta*dkap;
      if tn > 0 && kn > 0
        EVn = nsd_eval_all(cones, offs, xn);
        if EVn.ok, okstep = true; break; end
      end
      beta = beta / 2;
    end
    if ~okstep, break; end
    x = xn; y = y + beta*dy; s = s + beta*ds; tau = tn; kappa = kn; EV = EVn;
  end
end
res = struct('status', status, 'x', x, 'y', y, 's', s, 'tau', tau, 'kappa', kappa, ...
             'mu', muf(x, s, tau, kappa), 'iters', size(hist, 1), 'hist', hist, 'info', info, ...
             'time', toc(t0), 'EV', EV, 'offs', offs, 'last_cert', last_cert, 'cert_hist', {cert_hist}, 'last_fail', last_fail, 'accuracy', accuracy);
end

function p = prox(EV, offs, s, tau, kappa, mu)
psi = s + mu*EV.g;
p = sqrt(max(0, psi' * nsd_hsolve(EV, offs, psi)) + (tau*kappa - mu)^2) / mu;
end

function [dx, dy, dtau, ds, dkap] = direction(EV, offs, A, b, c, mu, tau, r1, r2, r3, rs, rk)
% Newton system of the homogeneous model, solved through the Schur complement (Lemma 3.1).
HiAt = nsd_hsolve(EV, offs, A');
Wm = HiAt / mu;
M = A * Wm; M = (M + M') / 2;
[Rm, p] = chol(M);
if p == 0
  Msol = @(v) Rm \ (Rm' \ v);
else
  Mi = pinv(M); Msol = @(v) Mi * v;
end
u = nsd_hsolve(EV, offs, r2 + rs) / mu;
v = nsd_hsolve(EV, offs, c) / mu;
pp = Msol(r1 - A*u);
qq = Msol(A*v + b);
w = b - Wm' * c;
dtau = (r3 + rk + c'*u - w'*pp) / (w'*qq + c'*v + mu / tau^2);
dy = pp + qq*dtau;
dx = u + Wm*dy - v*dtau;
ds = rs - mu * nsd_hmul(EV, offs, dx);
dkap = rk - (mu / tau^2) * dtau;
end
