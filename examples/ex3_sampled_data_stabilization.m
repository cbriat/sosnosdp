%% Example 3: sampled-data stabilization by clock-dependent LMIs, without SDP
%
% Plant x' = A x + B u with a zero-order hold, sampling period T. The augmented
% state (x, u) gives an impulsive system. We look for X(t), t in [0, T], and U with
%
%     dX/dt - Aa X - X Aa'  >= 0                       for all t in [0, T]
%
%     [ X(T)               , -(J0 X(T) + Ba U)' ]
%     [ -(J0 X(T) + Ba U)  ,  X(0) - eps I      ]  >= 0
%
% Then K = U X(T)^{-1} stabilizes the sampled loop. The first constraint is
% imposed directly on the interval (sosmatrixineq_interval). nsd_timevar lets us
% write the problem in physical time t, while SOSTOOLS stores it in the normalized
% variable s in [-1, 1], which keeps the monomial basis well conditioned.

% ---------------------------------------------------------------- data
A = [0.3 3; -3 0.3];          % unstable oscillator
B = [0; 1];

T = 2;                        % sampling period
epsl = 1e-3;                  % strictness margin

nx = size(A, 1);
nu = size(B, 2);
m = nx + nu;

% augmented matrices: state (x, u), u held constant between samples
Aa = [A B; zeros(nu, m)];
J0 = blkdiag(eye(nx), zeros(nu));
Ba = [zeros(nx, nu); eye(nu)];

% ---------------------------------------------------------------- degrees to try
for d = [6 8 32]

  % normalized clock variable s in [-1, 1] for t in [0, T]
  [s, tm] = nsd_timevar('s', [0 T]);

  prog = sosprogram(s);

  % decision variables: X(s) symmetric of degree d, U constant
  [prog, X] = sospolymatrixvar(prog, monomials(s, 0:d), [m m], 'symmetric');
  [prog, U] = sospolymatrixvar(prog, monomials(s, 0), [nu m]);

  % flow condition: dX/dt - Aa X - X Aa' >= 0 on the whole interval
  P = tm.ddt(X) - X*Aa' - (X*Aa')';
  prog = sosmatrixineq_interval(prog, P, s, tm.interval());

  % values of X at the ends of the sampling interval
  XT = tm.at(X, T);
  X0 = tm.at(X, 0);

  % jump condition (constant LMI)
  Off = -(XT*J0' + U'*Ba')';
  Psi = [XT,  Off';
         Off, X0 - epsl*eye(m)];
  prog = sosmatrixineq(prog, Psi, 'Mineq');

  % solve without SDP
  [prog, info] = sossolve_nosdp(prog);

  fprintf('d = %d: %s (%d iterations, %.2f s)\n', d, info.status, info.iters, info.time);

  % A candidate is returned also when the SOS certificate could not be established (status
  % 'uncertified_*', e.g. for badly conditioned plants). The closed-loop check below is
  % independent of the SOS certificate, so it can still validate the controller.
  certified = any(strcmp(info.status, {'feasible', 'optimal'}));
  candidate = ~isempty(strfind(info.status, 'uncertified'));

  if certified || candidate

    % controller gain K = U X(T)^{-1}
    XTv = double(sosgetsol(prog, XT));
    Uv = double(sosgetsol(prog, U));
    K = Uv / XTv;

    % independent check: spectral radius of the sampled closed loop
    Phi = expm(Aa * T);
    rho = max(abs(eig((J0 + Ba*K) * Phi)));

    fprintf('       K = %s\n', mat2str(K, 4));
    fprintf('       closed-loop spectral radius %.4f', rho);

    if ~certified && rho < 1
      fprintf('  (SOS certificate unavailable; stability verified by the spectral radius)');
    end
    fprintf('\n');
  end
end
