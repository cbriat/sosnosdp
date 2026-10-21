function out = nsd_hmul(EV, offs, v)
%NSD_HMUL  Multiply by the block-diagonal Hessian.
out = zeros(size(v));
for k = 1:numel(EV.H)
  sl = offs(k)+1:offs(k+1); out(sl) = EV.H{k} * v(sl);
end
end
