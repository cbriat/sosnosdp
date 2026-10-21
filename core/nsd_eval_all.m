function EV = nsd_eval_all(cones, offs, x, need_hess)
%NSD_EVAL_ALL  Evaluate every cone of a product at x. EV.ok, EV.g, EV.E{k}, EV.H{k}, EV.C{k} (chol of H).
if nargin < 4, need_hess = true; end
EV.ok = false; EV.g = zeros(size(x)); EV.E = cell(1, numel(cones)); EV.H = EV.E; EV.C = EV.E;
for k = 1:numel(cones)
  sl = offs(k)+1:offs(k+1);
  E = nsd_cone_eval(cones{k}, x(sl), need_hess);
  if ~E.ok, return; end
  EV.E{k} = E; EV.g(sl) = E.g;
  if need_hess
    [R, p] = chol(E.H);
    if p > 0
      Hs = E.H + 1e-12 * trace(E.H) / size(E.H, 1) * eye(size(E.H, 1));
      R = chol(Hs);
    end
    EV.H{k} = E.H; EV.C{k} = R;
  end
end
EV.ok = true;
end
