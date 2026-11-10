# TODO: scaling limits and ways around them

This file records where the current implementation stops scaling, how far each limit can be pushed
back, and whether certificates based on circuit polynomials (SONC/SAGE) would help. Line numbers
refer to the 0.1.0 code.

## 1. Current state

"Without SDP" holds in a precise but narrow sense:

* no SDP solver is called;
* SOS **constraints** are not lifted to Gram matrices; they are imposed through their values at N
  interpolation nodes, with the barrier `F(w) = -sum_j log det Lambda_j(w)`.

It does **not** mean "without Hessians" or "without PSD cones".

* The method is an interior-point method on non-symmetric cones. Each cone (N nodes, matrix size r)
  gets a **dense barrier Hessian** of order `dim = N r(r+1)/2`. It is rebuilt and Cholesky-factored
  at every line-search trial and corrector step (`core/nsd_cone_eval.m:17,35-37`,
  `core/nsd_eval_all.m:11`, `core/nsd_ipm.m:60,81`).
* Newton systems go through a **dense Schur complement** of order `n_f`, the number of free variables
  after elimination (`core/nsd_ipm.m:102-109`).
* Gram matrices of SOS constraints are never variables. They are reconstructed only to certify
  iterates (`core/nsd_gram.m`), at every iteration with tau > kappa in optimisation mode.
* Gram matrices of SOS **decision variables** (`sossosvar`, `sossosmatrixvar`, `sosposmatr`) are
  real variables. They enter as PSD cones, built as the special case N = 1, r = n of the interpolant
  cone (`bridge/sossolve_nosdp.m:58`).

### 1.1 Cost per iteration (dense, as implemented)

| step | time | memory |
|---|---|---|
| Hessian of one cone | O(r^4 N^2 + L r (N r)^2) | O(r^4 N^2); the vec temporary at `nsd_cone_eval.m:35` is ~4x the svec Hessian |
| Cholesky of that Hessian | O(dim^3) | O(dim^2) |
| Schur complement | O(sum_k dim_k^2 n_f + D n_f^2 + n_f^3) | O(D n_f + n_f^2) |

Here D is the total dimension of all cones.

For constraint cones, dim equals the number of coefficient constraints m. This is the same order as
the Schur complement of a Schur-complement SDP solver, which is the trade-off the manual describes.

### 1.2 Wall 1: Gram blocks of SOS decision variables

A Gram block of side n becomes a cone with dim = n(n+1)/2. Its entries are also free variables, so
they count in n_f. At n = 1800 (PIETOOLS-size blocks):

* dim = 1,620,900;
* the Hessian has dim^2 entries, about **21 TB**;
* the formation temporary has n^4 entries, about **84 TB**;
* one Cholesky costs about dim^3/3, about **1.4e18 flops**;
* the Schur complement is at least as large as the Hessian (n_f >= dim), about another 21 TB.

An SDP solver never forms this Hessian: it applies the PSD barrier Hessian through its X (x) X
structure in O(n^2) memory. Its Schur complement has the order of the number of coefficient
constraints, not of the number of Gram entries.

### 1.3 Wall 2: equality constraints

Equalities (`soseq`, and monomials outside the span of the Gram basis) are eliminated densely:

* `Ef = full(E)`, `y0 = pinv(Ef)*f` and `N = null(Ef)`, all based on a dense SVD
  (`bridge/sossolve_nosdp.m:114-119`);
* then `Lr = full(L*N)` (`bridge/sossolve_nosdp.m:123`);
* then a dense `qr(full(L), 0)` (`core/nsd_solve.m:14`).

E has one column per variable, including every svec entry of every Gram block, so the two walls
compound:

* 1e4 equality rows next to one 1800 Gram block make `Ef` about 130 GB, and `null(Ef)` about 21 TB;
* 1e5 rows by 1e5 variables is 80 GB for `Ef` alone, plus an SVD of order 1e15 flops.

### 1.4 Documentation that overstates the scaling

* [ ] `docs/manual.tex:185` claims "a memory footprint independent of the Gram blocks". This holds
  for SOS constraints only, not for Gram blocks of decision variables. Correct it.
* [ ] README "Limitations" and manual `sec:limits` say only "dense linear algebra". Add both walls,
  with the formulas of 1.1 and the 1800 example.

## 2. How far each limit can be circumvented

### 2.1 Gram blocks: largely avoidable, down to SDP parity

**Replace Gram variables by Gram-free interpolant variables.** The multiplier is represented by its
values at U nodes, with an interpolant SOS cone on them. This is what `sossosvar_nosdp` already does
for scalar multipliers.

* Variables and cone dimension drop from about (rL)^2/2 to U r(r+1)/2. Here L = C(n+d, n) is the
  size of the Gram basis and U = C(n+2d, n) the number of nodes.
* Scalar example: 3 variables, d = 20. L = 1771 and U = 12,341. The Hessian drops from about 20 TB
  to about 1.2 GB, and one Cholesky to about 6e11 flops.
* Matrix example: r = 10, 2 variables, d = 17. L = 171, so the Gram side is 1710, and U = 630. The
  cone dimension is 34,650 and the Hessian about 9.6 GB, against about 17 TB for the Gram block.
* What remains: the Hessian then has the order of the number of coefficient constraints, the same as
  the Schur complement of a Schur-complement SDP solver on the Gram formulation. This reaches
  **parity**, not an advantage, and dense matrices of that order are still large at PIETOOLS scale.

Tasks:

* [ ] `sossosmatrixvar_nosdp`: a Gram-free matrix multiplier. The core already supports matrix
  interpolant cones; only the bridge is missing.
* [ ] Warn in `sossolve_nosdp` when a Gram block of side above ~200 is present, and point to the
  Gram-free variants.

**Alternatives if Gram blocks must stay:**

* *Structured PSD cone oracle.* Apply `H[V] = X^-1 V X^-1` without forming H. This removes the n^4
  Hessian, but **not** the Schur complement over the n^2/2 Gram entries, as long as these remain free
  variables. Helping there would mean reformulating with Gram blocks on the primal side, with
  coefficient-matching equalities, as SDP solvers do. That is a substantial change to `nsd_ipm`.
* *DSOS/SDSOS* (Ahmadi and Majumdar): inner LP/SOCP approximations of the PSD constraint. Each 2x2
  block is a 3-dimensional cone, so the Hessian is block-diagonal with O(n^2) memory. But n_f still
  grows like n^2, so this only helps together with sparse Schur-complement linear algebra, and it is
  conservative.

### 2.2 Equalities: avoidable with sparse linear algebra

The dense SVD is an implementation choice, not a property of the method.

* [ ] **Sparse elimination.** Most `soseq` rows are sparse coefficient matchings. Eliminate one pivot
  variable per row with sparse LU/QR (SuiteSparse through MATLAB's sparse `qr`/`lu`), and keep the
  particular solution and the substitution sparse. Fill-in has to be monitored. A general null-space
  basis is dense, so do not compute one.
* [ ] **Fallback: keep the equalities inside the IPM** as a saddle-point system
  `[M E'; E 0]`. This avoids the SVD but does not shrink n_f.
* [ ] **Certificates under the fallback.** The current certificates rely on `E y = f` holding to
  rounding error, which elimination guarantees. Equalities kept inside the IPM would need a final
  sparse least-squares projection and a re-check. The re-check machinery already exists in
  `nsd_solve` (around lines 58-81).
* [ ] **Basis-independence reparametrisation.** Replace the dense `qr(full(L), 0)` in `nsd_solve.m:14`
  by a sparse rank-revealing QR, or by column scaling plus a rank test.

### 2.3 The floor: dense constraint Hessians and dense Schur complement

After 2.1 and 2.2, every constraint cone still has a dense Hessian of order m (its number of
coefficients), and the Schur complement of order n_f is dense. That matches dense Schur-complement
SDP solvers. Going below it requires structure:

* [ ] **Matrix-free Hessian.** The Hessian of an interpolant cone is applied through products with
  `P_j` and `Lambda_j^-1`, at about O(N L r^3) per product, without forming the dim x dim matrix.
  Hypatia is designed around Hessian *operators* of this kind. The Schur complement would then be
  solved iteratively (preconditioned CG). Ill-conditioning near the optimum is the known risk, and
  it is untested here.
* [ ] **Sparsity.** Correlative and term sparsity (e.g. TSSOS) split large Gram blocks into smaller
  ones in the SDP setting. The interpolant-cone analogue, splitting one constraint into several
  smaller cones, is an open item.
* Large-scale SDP is usually handled by first-order methods (SCS, CDCS, COSMO) or by Burer-Monteiro
  low-rank factorisation. For interpolant cones the required projection has no closed form, so
  neither transfers directly.

### 2.4 Cheap wins, independent of the above

* [ ] `nsd_gram.m:14,19`. `H \ Svec` and `H \ (Svec - Srec)` each re-factor the dense Hessian, up to
  four times per cone per certified iterate. Reuse `EV.C{k}`, the Cholesky factor already computed.
  This requires passing it through `nsd_check_values`.
* [ ] `nsd_cone_eval.m:35-37`. Build the svec Hessian directly instead of the vec temporary of order
  r^2 N, which saves about 4x memory.
* [ ] `nsd_ipm.m:60,81`. Every line-search trial builds and factors all Hessians. Test membership
  with `chol(Lambda_j)` alone first, and build Hessians only for candidate steps.
* [ ] In optimisation mode, certify every few iterations, or only once the gap is small, rather than
  at every iteration with tau > kappa.

## 3. Circuit polynomials (SONC / SAGE)

### 3.1 What they are

* A *circuit polynomial* has a support made of the vertices of a simplex (monomial squares with
  positive coefficients) plus at most one interior term. Its nonnegativity is decided in closed form
  by the "circuit number", with no Gram matrix (Iliman and de Wolff).
* Certificates are **sums of nonnegative circuit polynomials (SONC)**. They extend to compact
  semialgebraic sets through a Positivstellensatz with a converging hierarchy (Dressler, Iliman and
  de Wolff).
* They are computed by geometric or relative entropy programming. For polynomials, SONC coincides
  with the SAGE cone of Chandrasekaran and Shah (Murray, Chandrasekaran and Wierman).
* The SONC cone is also second-order-cone representable (Averkov; Wang and Magron), so SOCP solvers
  apply.
* Decompositions preserve sparsity: a SONC polynomial is a sum of nonnegative circuit polynomials
  supported on its own support (Wang; Murray et al.).
* Papp gives the duality theory and a column-generation algorithm for the circuits of an optimal
  bound, reporting sparse polynomials with up to 3000 monomials in minutes.

### 3.2 What they would solve here

* **No Gram blocks at all.** The size grows with the support and the number of circuits, not with
  (Gram side)^2. Wall 1 disappears for problems they can certify.
* **Native fit with the existing IPM.** Relative entropy cones are 3-dimensional exponential-type
  cones, so the Hessian is block-diagonal with tiny blocks. They are non-symmetric cones of exactly
  the kind the Skajaa-Ye method in `nsd_ipm` handles; alfonso and Hypatia do the same. A SAGE/SONC
  cone would be a new oracle in `core/`, without changing the solver.
* **Equalities (wall 2) are untouched.** That issue lives in our elimination code, not in the
  certificate.

### 3.3 Why they do not replace matrix SOS

* **Incomparable with SOS.** Neither cone contains the other (Iliman and de Wolff). SONC is strong
  for sparse polynomials of high degree with few terms, and can be weak or fail for dense supports,
  such as Lyapunov-type conditions with many cross terms.
* **Scalar only.** I found no established matrix-valued SONC/SAGE theory, and matrix SOS is the main
  purpose of this package. PIETOOLS constraints are matrix- or operator-valued.
* **Scalarising does not rescue it** (own argument, to check). One could scalarise a matrix
  constraint as `y' M(x) y`. But for quadratic forms in y, the only circuits are
  `{2e_i, 2e_j, e_i + e_j}`, so SONC in y amounts to scaled diagonal dominance, i.e. SDSOS-level
  conservatism in the matrix directions.

### 3.4 Verdict and tasks

SONC/SAGE is a useful **complementary scalar cone** for sparse, high-degree scalar constraints. It is
not a route to PIETOOLS-scale matrix programs.

* [ ] Optional prototype: a relative-entropy cone oracle in `core/` and `sosineq_sonc` in `bridge/`,
  tested against SOS on sparse examples, including the Motzkin polynomial (nonnegative, not SOS).

## 4. Suggested order

1. Documentation corrections (1.4).
2. Cheap wins (2.4).
3. `sossosmatrixvar_nosdp` and the large-Gram-block warning (2.1).
4. Sparse equality handling and sparse reparametrisation (2.2).
5. Matrix-free Hessians with iterative Schur solves (2.3), as a research item.
6. SONC/SAGE cone (3.4), optional.

## References

* D. Papp, S. Yildiz. Sum-of-squares optimization without semidefinite programming. SIAM J. Optim.
  29(1):822-851, 2019.
* L. Kapelevich, C. Coey, J. P. Vielma. Sum of squares generalizations for conic sets. Math. Program.,
  2023. doi:10.1007/s10107-022-01831-6.
* C. Coey, L. Kapelevich, J. P. Vielma. Performance enhancements for a generic conic interior point
  algorithm. Math. Program. Comput., 2023. [arXiv:2005.01136](https://arxiv.org/abs/2005.01136).
* S. Iliman, T. de Wolff. Amoebas, nonnegative polynomials and sums of squares supported on circuits.
  Res. Math. Sci., 2016. [arXiv:1402.0462](https://arxiv.org/abs/1402.0462).
* M. Dressler, S. Iliman, T. de Wolff. A Positivstellensatz for sums of nonnegative circuit
  polynomials. [arXiv:1607.06010](https://arxiv.org/abs/1607.06010).
* V. Chandrasekaran, P. Shah. Relative entropy relaxations for signomial optimization. SIAM J.
  Optim., 2016. [arXiv:1409.7640](https://arxiv.org/abs/1409.7640).
* R. Murray, V. Chandrasekaran, A. Wierman. Newton polytopes and relative entropy optimization.
  Found. Comput. Math., 2021. [arXiv:1810.01614](https://arxiv.org/abs/1810.01614).
* R. Murray, V. Chandrasekaran, A. Wierman. Signomial and polynomial optimization via relative
  entropy and partial dualization. Math. Program. Comput.
  [arXiv:1907.00814](https://arxiv.org/abs/1907.00814).
* J. Wang. Nonnegative polynomials and circuit polynomials. SIAM J. Appl. Algebra Geom. 6, 2022.
  doi:10.1137/20M1313969. [arXiv:1804.09455](https://arxiv.org/abs/1804.09455).
* J. Wang, V. Magron. A second order cone characterization for sums of nonnegative circuits.
  [arXiv:1906.06179](https://arxiv.org/abs/1906.06179). Follow-up:
  [arXiv:2012.07903](https://arxiv.org/abs/2012.07903).
* D. Papp. Duality of sum of nonnegative circuit polynomials and optimal SONC bounds.
  [arXiv:1912.04718](https://arxiv.org/abs/1912.04718).
* M. Dressler, H. Naumann, T. Theobald. The dual cone of sums of non-negative circuit polynomials.
  [arXiv:1809.07648](https://arxiv.org/abs/1809.07648).
* A. A. Ahmadi, A. Majumdar. DSOS and SDSOS optimization: more tractable alternatives to sum of
  squares and semidefinite optimization. SIAM J. Appl. Algebra Geom. 3(2):193-230, 2019.
  [arXiv:1706.02586](https://arxiv.org/abs/1706.02586).
* J. Wang, V. Magron, J.-B. Lasserre. TSSOS: a moment-SOS hierarchy that exploits term sparsity.
  SIAM J. Optim. 31(1):30-58, 2021. [arXiv:1912.08899](https://arxiv.org/abs/1912.08899).
