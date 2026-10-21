function [prog, info] = sossolve_nosdp(prog, opts)
%SOSSOLVE_NOSDP  Solve a SOSTOOLS program directly in matrix weighted SOS cones, without SDP.
%
%   [prog, info] = sossolve_nosdp(prog)  replaces  prog = sossolve(prog).
%   Afterwards sosgetsol(prog, V) works as usual.
%
%   Supported constraints
%     sosineq(prog, p)                 scalar SOS (global)
%     sosmatrixineq(prog, M)           matrix SOS, both 'quadraticMineq' (default) and 'Mineq'
%     sosineq(prog, p, [a inf])        half-line constraints (SOSTOOLS 'posint')
%     soseq(prog, p)                   equality constraints (eliminated)
%     sosineq_interval / sosmatrixineq_interval   finite intervals, natively weighted
%   Helpers: sossosvar_nosdp (SOS polynomial without Gram matrix), nsd_timevar (normalized variable)
%   Supported decision variables: sosdecvar, sospolyvar, sospolymatrixvar (free),
%   sossosvar_nosdp (SOS polynomial without Gram matrix), and sossosvar, sossosmatrixvar,
%   sosposmatr (their Gram matrices, which SOSTOOLS creates, enter as PSD cones of the same
%   method; no SDP solver is called).
%
%   opts (all optional): max_iter, tol, verbose, box ([lo hi] for node selection, default [-1 1]),
%     certificate ('auto' default and 'gram': explicit Gram matrices decide, used only for
%                  verification; 'dikin': Gram-free, faster, less robust near the boundary),
%     rank_tol (1e-14), depth_tol (1e-8): depth a Dikin certificate must prove
%   info: status ('feasible','optimal','infeasible',...), obj, iters, time, cones, gram, farkas_radius
if nargin < 2, opts = struct(); end
if ~isempty(prog.solinfo.x), error('The SOS program is already solved.'); end
box = [-1 1]; if isfield(opts, 'box'), box = opts.box; opts = rmfield(opts, 'box'); end
nU = prog.var.idx{end} - 1;
nvt = numel(prog.vartable);

% ---------------------------------------------------------------- decision variables
% u (SOSTOOLS vector) = Tmap * y; free entries map one to one, Gram blocks of 'sos' variables
% are parametrized by svec, so symmetric pairs and SOSTOOLS' split of cross terms are handled.
free = 1:prog.var.idx{1}-1;
psd = zeros(0, 2);                                   % [start, n]
for v = 1:prog.var.num
  i0 = prog.var.idx{v}; len = prog.var.idx{v+1} - i0;
  if strcmp(prog.var.type{v}, 'sos')
    n = round(sqrt(len)); assert(n*n == len, 'unexpected Gram block size');
    psd(end+1, :) = [i0, n]; %#ok<AGROW>
  else
    free = [free, i0:i0+len-1]; %#ok<AGROW>
  end
end
nf = numel(free);
ny = nf + sum(psd(:, 2) .* (psd(:, 2) + 1) / 2);
Ti = free(:); Tj = (1:nf)'; Tv = ones(nf, 1); col = nf;
for k = 1:size(psd, 1)
  n = psd(k, 2); S = nsd_svec_basis(n); [ii, jj, vv] = find(S);
  Ti = [Ti; psd(k, 1) - 1 + ii]; Tj = [Tj; col + jj]; Tv = [Tv; vv]; col = col + size(S, 2);
end
Tmap = sparse(Ti, Tj, Tv, nU, ny);

cones = {}; Lc = {}; cc = {}; E = sparse(0, ny); f = zeros(0, 1); desc = {};
% Gram blocks of SOS decision variables: PSD cones on their svec
col = nf;
for k = 1:size(psd, 1)
  n = psd(k, 2); t = n*(n+1)/2;
  cones{end+1} = nsd_cone(n, {1}, sprintf('Gram of SOS variable %d', k)); %#ok<AGROW>
  Lc{end+1} = sparse(1:t, col + (1:t), 1, t, ny); cc{end+1} = zeros(t, 1); %#ok<AGROW>
  desc{end+1} = cones{end}.name; col = col + t; %#ok<AGROW>
end

% ---------------------------------------------------------------- constraints
for e = 1:prog.expr.num
  At = prog.expr.At{e}; b0 = full(prog.expr.b{e}(:)); Z = full(prog.expr.Z{e});
  if size(At, 1) < nU, At = [At; sparse(nU - size(At, 1), size(At, 2))]; end %#ok<AGROW>
  C = -(At' * Tmap);                                 % coefficients = b0 + C*y
  nZ = size(Z, 1); dimp = round(sqrt(numel(b0) / nZ));
  type = prog.expr.type{e};
  if strcmp(type, 'eq')
    E = [E; C]; f = [f; -b0]; continue; %#ok<AGROW>
  end
  if size(Z, 2) > nvt && strcmp(type, 'sparsemultipartite')
    [Z, b0, C, dimp] = multipartite_to_matrix(Z, b0, C, nvt);
    nZ = size(Z, 1);
  end
  % active monomials (some coefficient can be nonzero)
  act = any(reshape(b0 ~= 0, nZ, []), 2) | any(reshape(any(C ~= 0, 2), nZ, []), 2);
  if ~any(act), continue; end
  ka = find(act);
  idx = bsxfun(@plus, ka, (0:dimp^2-1) * nZ); idx = idx(:);
  Za = Z(ka, :); b0 = b0(idx); C = C(idx, :); nZ = numel(ka);
  switch type
    case {'ineq', 'sparse', 'sparsemultipartite'}
      gb = [];
      if isfield(prog, 'nosdp') && isfield(prog.nosdp, 'gram_basis') && numel(prog.nosdp.gram_basis) >= e
        gb = prog.nosdp.gram_basis{e};
      end
      [Kc, cv, Lv, Ee, fe, dsc] = global_cone(Za, dimp, b0, C, box, gb);
    case 'posint'
      iv = find(any(Za, 1)); assert(numel(iv) <= 1, 'posint: univariate polynomial expected');
      if isempty(iv), iv = 1; end
      [Kc, cv, Lv, dsc] = interval_cone(Za, dimp, b0, C, iv, [0 inf]); Ee = []; fe = [];
    case 'nosdp_interval'
      itv = prog.nosdp.interval{e};
      iv = find(strcmp(prog.vartable, itv.var));
      assert(~isempty(iv), 'interval variable not found in the program variables');
      others = setdiff(find(any(Za, 1)), iv);
      assert(isempty(others), 'interval constraint must depend on its variable only');
      [Kc, cv, Lv, dsc] = interval_cone(Za, dimp, b0, C, iv, itv.ab); Ee = []; fe = [];
    otherwise
      error('sossolve_nosdp: constraint type ''%s'' is not supported', type);
  end
  cones{end+1} = Kc; Lc{end+1} = Lv; cc{end+1} = cv; desc{end+1} = dsc; %#ok<AGROW>
  if ~isempty(Ee), E = [E; Ee]; f = [f; fe]; end %#ok<AGROW>
end
L = vertcat(Lc{:}); c = vertcat(cc{:});

% ---------------------------------------------------------------- equalities and objective
obj = zeros(nU, 1);
if ~isempty(prog.objective), obj(1:numel(prog.objective)) = full(prog.objective); end
bobj = -(Tmap' * obj);                               % SOSTOOLS minimizes obj'u
if size(E, 1) > 0
  Ef = full(E);
  y0 = pinv(Ef) * f;
  if norm(Ef * y0 - f) > 1e-8 * max(1, norm(f))
    info = struct('status', 'infeasible_equalities', 'iters', 0); prog.solinfo.info = info; return;
  end
  N = null(Ef);
else
  y0 = zeros(ny, 1); N = speye(ny);
end
cr = c + L * y0; Lr = full(L * N); br = N' * bobj;

% ---------------------------------------------------------------- solve
out = nsd_solve(cones, Lr, cr, br, opts);
info = struct('status', out.status, 'iters', out.iters, 'time', out.time, 'cones', {desc}, ...
              'gram', {out.gram}, 'n_free', size(Lr, 2), 'cone_dim', size(Lr, 1), ...
              'barrier_parameter', sum(cellfun(@(K) K.nu, cones)), 'accuracy', out.accuracy);
if isfield(out, 'farkas_radius'), info.farkas_radius = out.farkas_radius; end
if isfield(out, 'obj_gap'), info.obj_gap_to_converged = out.obj_gap; end
info.rank = out.rank;
info.cert_method = out.cert_method; info.dikin_ratio = out.dikin_ratio;
if isfield(out, 'depth'), info.depth = out.depth; end
if isfield(out, 'cert_failure'), info.cert_failure = out.cert_failure; end
info.last_failed_certificate = out.ipm.last_fail;
info.near_dependent = out.near_dependent;
if out.near_dependent > 0
  warning(['sossolve_nosdp: %d directions of the decision variables are numerically dependent ', ...
           '(rank %d of %d). This usually comes from a monomial basis on a wide interval; ', ...
           'declare the independent variable with nsd_timevar (normalized to [-1, 1]).'], ...
          out.near_dependent, out.rank, size(Lr, 2));
end
if isfield(out, 'coef_consistency'), info.coef_consistency = out.coef_consistency; end
if ~isempty(out.z)
  y = y0 + N * out.z;
  u = Tmap * y;
  prog.solinfo.RRx = u; prog.solinfo.x = u;
  info.obj = obj' * u;
end
prog.solinfo.info = info;
end

% ======================================================================== helpers
function [Kc, cv, Lv, Ee, fe, dsc] = global_cone(Za, dimp, b0, C, box, gb)
% Matrix SOS on R^n: basis monomials (given, or from the Newton box), nodes unisolvent for products.
if nargin >= 6 && ~isempty(gb)
  vars = find(any(Za, 1) | any(gb, 1)); if isempty(vars), vars = 1; end
  Zv = Za(:, vars); Zh = gb(:, vars);
else
  vars = find(any(Za, 1)); if isempty(vars), vars = 1; end
  Zv = Za(:, vars); deg = sum(Zv, 2);
  Zh = nsd_monomials(numel(vars), ceil(min(deg) / 2), floor(max(deg) / 2));
  lo = ceil(min(Zv, [], 1) / 2); hi = floor(max(Zv, [], 1) / 2);
  Zh = Zh(all(bsxfun(@ge, Zh, lo), 2) & all(bsxfun(@le, Zh, hi), 2), :);
end
assert(~isempty(Zh), 'sossolve_nosdp: no Gram monomials (odd or empty polynomial)');
nh = size(Zh, 1);
[I, J] = ndgrid(1:nh, 1:nh);
ZZ = unique(Zh(I(:), :) + Zh(J(:), :), 'rows');
% monomials outside span(ZZ) must vanish
[inZZ] = ismember(Zv, ZZ, 'rows');
Ee = []; fe = [];
if any(~inZZ)
  bad = find(~inZZ); nZ = size(Zv, 1);
  idx = bsxfun(@plus, bad, (0:dimp^2-1) * nZ); idx = idx(:);
  Ee = C(idx, :); fe = -b0(idx);
end
X = nsd_nodes(ZZ, box(1) * ones(1, numel(vars)), box(2) * ones(1, numel(vars)));
P = nsd_orth(nsd_moneval(Zh, X));
Kc = nsd_cone(dimp, {P}, sprintf('SOS, %d vars, matrix size %d, %d nodes', numel(vars), dimp, size(X, 1)));
[cv, Lv] = nsd_values_map(Zv, dimp, b0, C, X, 1:numel(vars), Kc.S);
dsc = Kc.name;
end

function [Kc, cv, Lv, dsc] = interval_cone(Za, dimp, b0, C, iv, ab)
% Univariate (matrix) weighted SOS on [a,b] (weights 1, (x-a)(b-x)) or [a,inf) (weights 1, x-a).
D = max(Za(:, iv));
a = ab(1); b = ab(2);
if isinf(b)
  n0 = floor(D / 2); n1 = floor((D - 1) / 2); Nn = D + 1;
  xs = a + (1 - cos(pi * (0:Nn-1)' / max(Nn - 1, 1)));          % nodes in [a, a+2]
  if Nn == 1, xs = a + 1; end
  wgt = xs - a; xm = xs - a - 1;                                 % Chebyshev variable on [a, a+2]
else
  n0 = ceil(D / 2); n1 = n0 - 1; Nn = 2 * n0 + 1;
  xm = cos(pi * (0:Nn-1)' / max(Nn - 1, 1)); if Nn == 1, xm = 0; end
  xs = (a + b) / 2 + (b - a) / 2 * xm;
  wgt = (xs - a) .* (b - xs);
end
Ps = {nsd_orth(cheb(xm, n0))};
if n1 >= 0
  Ps{2} = nsd_orth(bsxfun(@times, sqrt(max(wgt, 0)), cheb(xm, n1)));
end
Kc = nsd_cone(dimp, Ps, sprintf('interval [%g, %g], degree %d, matrix size %d', a, b, D, dimp));
[cv, Lv] = nsd_values_map(Za(:, iv), dimp, b0, C, xs, 1, Kc.S);
dsc = Kc.name;
end

function T = cheb(x, n)
% Chebyshev polynomials T_0..T_n at x (column vector).
T = zeros(numel(x), n + 1); T(:, 1) = 1;
if n >= 1, T(:, 2) = x; end
for k = 3:n+1, T(:, k) = 2 * x .* T(:, k-1) - T(:, k-2); end
end

function [Zx, b1, C1, dimp] = multipartite_to_matrix(Z, b0, C, nvt)
% y' M(x) y with auxiliary variables y (columns > nvt) -> matrix polynomial M(x).
Zxa = Z(:, 1:nvt); Zy = Z(:, nvt+1:end);
assert(all(sum(Zy, 2) == 2 | (b0 == 0 & ~any(C, 2))), 'multipartite expression is not quadratic in y');
used = find(any(Zy, 1)); dimp = numel(used); Zy = Zy(:, used);
[Zx, ~, kx] = unique(Zxa, 'rows'); nx = size(Zx, 1);
b1 = zeros(nx * dimp^2, 1); C1 = sparse(nx * dimp^2, size(C, 2));
for k = 1:size(Z, 1)
  pq = find(Zy(k, :));
  if numel(pq) == 1, p = pq; q = pq; w = 1; else, p = pq(1); q = pq(2); w = 0.5; end
  for ent = unique([(q-1)*dimp + p, (p-1)*dimp + q])
    row = (ent - 1) * nx + kx(k);
    b1(row) = b1(row) + w * b0(k); C1(row, :) = C1(row, :) + w * C(k, :);
  end
end
end
