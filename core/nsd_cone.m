function K = nsd_cone(r, Ps, name)
%NSD_CONE  Matrix weighted SOS cone in an interpolant basis (dual barrier -sum log det Lambda_j).
%   K = nsd_cone(r, Ps, name): r = matrix size, Ps = cell of N-by-L_j matrices whose rows are
%   sqrt(g_j(xi_i)) * b_j(xi_i)' (weights folded in), N = number of nodes.
%   A point of the cone is stored node by node: [svec(S_1); ...; svec(S_N)].
%   The PSD cone of size r is nsd_cone(r, {1}).
if nargin < 3, name = ''; end
if ~iscell(Ps), Ps = {Ps}; end
K.r = r; K.t = r*(r+1)/2; K.N = size(Ps{1}, 1); K.Ps = Ps; K.name = name;
K.dim = K.N * K.t;
K.L = cellfun(@(P) size(P, 2), Ps);
K.nu = r * sum(K.L);
K.S = nsd_svec_basis(r);
K.Sbig = kron(speye(K.N), sparse(K.S));
K.Pk = cell(size(Ps));
for j = 1:numel(Ps)
  K.Pk{j} = kron(Ps{j}, eye(r));          % (N r) x (L_j r)
end
% reference interior element E = A(I) has node values (sum_j ||p_ji||^2) I; its largest entry
Enode = zeros(K.N, 1);
for j = 1:numel(Ps), Enode = Enode + sum(Ps{j}.^2, 2); end
K.Enorm = max(Enode);
% initial point: identity at every node
K.x0 = repmat(K.S' * reshape(eye(r), [], 1), K.N, 1);
end
