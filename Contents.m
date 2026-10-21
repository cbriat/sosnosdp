% SOSNOSDP  Matrix sum-of-squares programs in SOSTOOLS, solved without semidefinite programming.
% Version 0.1.0
%
% Setup
%   setup_nosdp             - add SOSNOSDP to the path and check the SOSTOOLS installation
%
% Solving SOSTOOLS programs
%   sossolve_nosdp          - replacement for sossolve; sosgetsol works afterwards
%   sosineq_interval        - scalar polynomial inequality on an interval [a, b]
%   sosmatrixineq_interval  - matrix polynomial inequality on an interval [a, b]
%   sossosvar_nosdp         - SOS decision polynomial without a Gram matrix
%   nsd_timevar             - normalized variable on [-1, 1] for a physical interval
%
% Low-level solver (core/)
%   nsd_cone                - matrix weighted SOS cone in an interpolant basis
%   nsd_cone_eval           - barrier, gradient and Hessian of the dual cone
%   nsd_solve               - max b'z s.t. c + L z in K, with certificates
%   nsd_ipm                 - homogeneous predictor-corrector interior-point method
%
% See README.md and docs/manual.pdf.
