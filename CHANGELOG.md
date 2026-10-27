# Changelog

## 0.1.0 (first public version)

* `sossolve_nosdp`: solves SOSTOOLS 4 programs (`pvar` mode) in interpolant matrix weighted SOS cones
  with a homogeneous interior-point method; no SDP solver and no Gram lifting of the constraints.
* Supported: `sosineq`, `sosmatrixineq` (both forms), half-line constraints, `soseq`, `sosdecvar`,
  `sospolyvar`, `sospolymatrixvar`, `sossosvar`, `sossosmatrixvar`, `sosposmatr`, `sossetobj`, `sosgetsol`.
* New functions: `sosineq_interval`, `sosmatrixineq_interval`, `sossosvar_nosdp`, `nsd_timevar`.
* Certificates for the returned coefficients: explicit Gram matrices with refinement and rounding
  checks (default), a Gram-free Dikin test as an option, Farkas vectors with a radius for infeasible
  problems.
* Five examples, 21 tests, Octave support (with documented edits to SOSTOOLS), CI workflow.

Fixed before release:

* a near-boundary point could be accepted as a certificate (regression test T12);
* the Dikin test accepted, with computed depth 1e-8, a two-variable matrix polynomial that is not
  positive semidefinite far from the nodes (T15). The Gram test now decides in the default mode, with
  iterative refinement of the Gram identity and a residual required to be small relative to the
  proved margin; the last ten certified iterates are kept as fallbacks.
* the manual compared the cost per iteration with an SDP accounting that applies only to full-KKT
  solvers; it now gives the counts for Schur-complement and full-KKT solvers separately.
* the Gram certificate failed with `MATLAB:ndFullOnly` in `nsd_gram` when a cone had dimension one
  (the Gram block of a single-monomial `sossosvar` or `sossosmatrixvar`, or a scalar SOS constraint
  interpolated at a single node, e.g. `sosineq(prog, 1 - a)`): a sparse 1x1 product was indexed as a
  4-D array. The Gram matrices and that product are now full; the `'dikin'` certificate was not
  affected (regression test T18).
