function out = nsd_hsolve(EV, offs, V)
%NSD_HSOLVE  Apply the inverse of the block-diagonal Hessian to the columns of V.
out = zeros(size(V));
for k = 1:numel(EV.C)
  sl = offs(k)+1:offs(k+1); R = EV.C{k};
  out(sl, :) = R \ (R' \ V(sl, :));
end
end
