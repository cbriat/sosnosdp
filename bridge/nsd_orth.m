function [Q, R] = nsd_orth(P)
%NSD_ORTH  Orthonormalize the columns of P (Papp-Yildiz, Sec. 4.1): P = Q R.
[Q, R] = qr(P, 0);
end
