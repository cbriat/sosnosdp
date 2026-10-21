# Changelog

## 0.1.0 (first public version)

* `sossolve_nosdp`: solves SOSTOOLS 4 programs (`pvar` mode) in interpolant matrix weighted SOS cones
  with a homogeneous interior-point method; no SDP solver and no Gram lifting of the constraints.
* Supported: `sosineq`, `sosmatrixineq` (both forms), half-line constraints, `soseq`, `sosdecvar`,
  `sospolyvar`, `sospolymatrixvar`, `sossosvar`, `sossosmatrixvar`, `sosposmatr`, `sossetobj`, `sosgetsol`.
* New functions: `sosineq_interval`, `sosmatrixineq_interval`, `sossosvar_nosdp`, `nsd_timevar`.
* Certificates for the returned coefficients: Dikin test with a scale-free depth requirement, explicit
  Gram matrices as verification, Farkas vectors with a radius for infeasible problems.
* Five examples, sixteen tests, Octave support (with documented edits to SOSTOOLS), CI workflow.

Fixed before release: a near-boundary point could be accepted as a certificate (regression test T12).
