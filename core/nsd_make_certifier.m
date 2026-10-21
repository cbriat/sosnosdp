function cert = nsd_make_certifier(cones, L, c, mode, depth_tol)
%NSD_MAKE_CERTIFIER  Certificate test for s(z) = c + L z in K at an iterate (see NSD_CERTIFY_POINT).
if nargin < 4, mode = 'auto'; end
if nargin < 5, depth_tol = 1e-8; end
offs = [0 cumsum(cellfun(@(K) K.dim, cones))];
cert = @(x, y, tau, kappa, mu, EV) nsd_certify_point(cones, offs, L, c, mu / tau, EV, mode, depth_tol, x);
end
