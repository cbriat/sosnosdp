function [ok, info] = nsd_certify_point(cones, offs, L, c, delta, EV, mode, depth_tol, x)
%NSD_CERTIFY_POINT  Certificate that c + L z* lies in the interior of the cone K, at an iterate.
%   z* = argmin ||c + L z + delta g||_{H^-1} (weighted projection, eq. (5.1) of the report), then
%     'dikin' : ||S_k + delta g_k||_{H_k^-1} < delta for every cone k, i.e. S_k lies in the Dikin
%               ellipsoid of the conjugate barrier (Theorem 5.2's sufficient condition). Uses only
%               the barrier Hessian; no Gram matrix is formed.
%     'gram'  : explicit Gram matrices Phi_j = Lambda_j(x)^-1 Lambda_j(w) Lambda_j(x)^-1 (w = H^-1 S),
%               accepted if all are positive definite.
%     'auto'  : as 'gram' (default); the Dikin ratio is computed and reported only.
%   Acceptance rules are in NSD_CHECK_VALUES (the Dikin test must prove depth >= depth_tol).
if nargin < 7 || isempty(mode), mode = 'auto'; end
if nargin < 8 || isempty(depth_tol), depth_tol = 1e-8; end
if nargin < 9, x = []; end
ok = false; info = struct();
if ~(delta > 0 && isfinite(delta)), return; end
HiL = nsd_hsolve(EV, offs, L);
Hir = nsd_hsolve(EV, offs, c + delta * EV.g);
Mw = L' * HiL; Mw = (Mw + Mw') / 2;
z = -(Mw \ (L' * Hir));
S = c + L*z;
[ok, info] = nsd_check_values(cones, offs, S, delta, EV, mode, depth_tol);
info.z = z; info.x = x;
end
