function [msq, rho_periodic] = second_moment(sys, Kd, x0, thetas)
%SECOND_MOMENT  Exact second moment of the stochastic sampled-data system (no Monte Carlo).
%
%   Between samples, Sigma = E[x x'] obeys Sigma' = Abar Sigma + Sigma Abar' + E1 Sigma E1'
%   + E2 Sigma E2'; at a sample, Sigma+ = Jcl Sigma Jcl' with Jcl = J0 + Bbar Kd. In vectorized
%   form (Theorem 5 (ii) of the paper) this is linear, so the second moment after any given
%   sequence of sampling intervals is computed exactly.
%
%   msq(k)         : E||x(t_k^+)||^2 after the k-th sample, for the intervals in thetas
%   rho_periodic(i): spectral radius of the periodic second-moment map for interval thetas(i)

m = size(sys.Abar, 1);
Im = eye(m);

% generator of vec(Sigma) between samples
Acal = kron(Im, sys.Abar) + kron(sys.Abar, Im) + kron(sys.E1, sys.E1) + kron(sys.E2, sys.E2);

% jump map of vec(Sigma)
Jcl = sys.J0 + sys.Bbar * Kd;
Jcal = kron(Jcl, Jcl);

v = reshape(x0 * x0', [], 1);
msq = zeros(numel(thetas), 1);
rho_periodic = zeros(numel(thetas), 1);

for k = 1:numel(thetas)
  Mk = Jcal * expm(Acal * thetas(k));
  v = Mk * v;
  msq(k) = trace(reshape(v, m, m));
  rho_periodic(k) = max(abs(eig(Mk)));
end
end
