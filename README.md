# SOSNOSDP

**Matrix sum-of-squares programs in SOSTOOLS, solved without semidefinite programming.**

SOSNOSDP is a MATLAB package (also running in GNU Octave) that solves SOS programs written with
[SOSTOOLS](https://github.com/oxfordcontrol/SOSTOOLS) without any SDP solver. Scalar and matrix
polynomial constraints are represented by their values at interpolation nodes, and the program is solved
by a non-symmetric interior-point method that only evaluates the barrier of the dual cone. No Gram
matrices are introduced for the constraints, and every returned solution is certified.

Replace

```matlab
prog = sossolve(prog);
```

by

```matlab
[prog, info] = sossolve_nosdp(prog);
```

and use `sosgetsol` as usual.

## Built on

The method is not new. This package makes it available to SOSTOOLS users.

* **D. Papp and S. Yıldız**, *Sum-of-squares optimization without semidefinite programming*,
  SIAM J. Optim. 29 (2019) 822–851: SOS optimization in the interpolant basis with a non-symmetric
  interior-point method, the dual barrier and the recovery of SOS decompositions. Their MATLAB solver is
  [alfonso](https://github.com/dpapp-github/alfonso) (INFORMS J. Comput. 34 (2022) 11–19).
* **L. Kapelevich, C. Coey and J. P. Vielma**, *Sum of squares generalizations for conic sets*,
  Math. Program., [doi:10.1007/s10107-022-01831-6](https://doi.org/10.1007/s10107-022-01831-6):
  the cone of SOS matrices in the interpolant basis, its barrier and oracles, building on
  D. Papp and F. Alizadeh, SIAM J. Optim. 23 (2013) 1398–1423.
* **C. Coey, L. Kapelevich and J. P. Vielma**, *Solving natural conic formulations with Hypatia.jl*,
  INFORMS J. Comput. 34 (2022) 2686–2699: the Julia solver [Hypatia](https://github.com/jump-dev/Hypatia.jl),
  which implements these cones (`WSOSInterpPosSemidefTriCone`).
* **A. Skajaa and Y. Ye**, Math. Program. 150 (2015) 391–422: the homogeneous interior-point method.
* **SOSTOOLS** (A. Papachristodoulou, J. Anderson, G. Valmorbida, S. Prajna, P. Seiler, P. A. Parrilo,
  M. M. Peet, D. Jagt): the modelling front end.

What SOSNOSDP adds:

* a compact MATLAB implementation of the matrix weighted SOS cone and of the interior-point method,
* the SOSTOOLS bridge `sossolve_nosdp`,
* interval constraints for scalar and matrix polynomials, without multipliers,
* a helper for normalized time variables,
* certificates for the coefficients actually returned.

In Julia, [SumOfSquares.jl](https://github.com/jump-dev/SumOfSquares.jl) with Hypatia offers SDP-free
scalar SOS constraints. See the manual for a comparison.

## Installation

Requirements: MATLAB R2016b or later (or GNU Octave 8) and SOSTOOLS 4, with programs written using
`pvar`. No SDP solver is needed.

```matlab
addpath(genpath('/path/to/SOSTOOLS'));
addpath('/path/to/sosnosdp');
setup_nosdp
```

For GNU Octave, SOSTOOLS needs three one-line edits, applied by `octave/patch_sostools.py`
(see [octave/README_OCTAVE.md](octave/README_OCTAVE.md)).

## Quick start

```matlab
pvar x a

M = [1 + x^2, a*x;
     a*x,     2 + x^2];

prog = sosprogram(x, a);
prog = sosmatrixineq(prog, M);      % M(x) positive semidefinite for all x
prog = sossetobj(prog, -a);         % maximize a

[prog, info] = sossolve_nosdp(prog);

double(sosgetsol(prog, a))          % 2.4142..., the exact value is 1 + sqrt(2)
```

## What is supported

| SOSTOOLS construct | Treatment |
|---|---|
| `sosineq(prog, p)` | interpolant SOS cone |
| `sosmatrixineq(prog, M)` (default and `'Mineq'`) | matrix SOS cone |
| `sosineq(prog, p, [a inf])` | weighted cone on a half-line |
| `soseq` | linear equalities, eliminated |
| `sosdecvar`, `sospolyvar`, `sospolymatrixvar` | free decision variables |
| `sossosvar`, `sossosmatrixvar`, `sosposmatr` | Gram blocks created by SOSTOOLS, as PSD cones of the same method |
| `sossetobj`, `sosgetsol` | as in SOSTOOLS |

New functions:

| Function | Purpose |
|---|---|
| `sosineq_interval(prog, p, t, [a b])` | `p(t) >= 0` on `[a, b]`, no multiplier |
| `sosmatrixineq_interval(prog, M, t, [a b])` | `M(t)` positive semidefinite on `[a, b]`, no multiplier |
| `sossosvar_nosdp(prog, Z)` | SOS decision polynomial without a Gram matrix |
| `[s, tm] = nsd_timevar('s', [a b])` | normalized variable on `[-1, 1]`, with `tm.ddt`, `tm.at`, `tm.interval` |

## Examples

| File | Content |
|---|---|
| `examples/ex1_polynomial_minimum.m` | global and interval minimum of a polynomial |
| `examples/ex2_matrix_sos.m` | matrix SOS |
| `examples/ex3_sampled_data_stabilization.m` | sampled-data stabilization by clock-dependent LMIs |
| `examples/ex4_aperiodic_stability.m` | largest sampling interval under aperiodic sampling, by bisection |
| `examples/ex5_stochastic_sampled_data.m` | mean-square stabilization of a stochastic sampled-data system (C. Briat, Automatica 74 (2016), Theorem 15) |

## Documentation

The manual [docs/manual.pdf](docs/manual.pdf) covers:

* how the method works: the interpolant representation, the matrix WSOS cone and its barrier, the
  interior-point method, the certificates, and the translation of SOSTOOLS programs,
* all functions, options and statuses,
* practical advice on conditioning,
* a walk through the examples.

Two pieces of advice from it:

* **Declare interval variables with `nsd_timevar`.** SOSTOOLS stores polynomials in monomial
  coefficients, and on wide intervals high-degree monomials become numerically dependent.
* **Prefer `sosmatrixineq_interval` to SOS multipliers** for interval constraints.

## Tests

`tests/run_tests.m` runs 16 tests with known answers. They pass under GNU Octave 8.4 with SOSTOOLS 4,
and run automatically on every push (`.github/workflows/octave-tests.yml`). Comparisons with SeDuMi
through `sossolve` are in `tests/compare_sedumi.m` and `tests/scaling_sedumi.m`.

The package has not yet been run in MATLAB by its authors. Reports are welcome.

## Limitations

* Only SOSTOOLS programs written with `pvar`, not the Symbolic Math Toolbox mode.
* Sparsity options of SOSTOOLS are handled as dense constraints.
* Global SOS constraints in several variables are much less tested than univariate and interval ones.
* Dense linear algebra, and certificates evaluated in floating point.

## How to cite

Please cite this package (see [CITATION.cff](CITATION.cff)) together with Papp and Yıldız (2019),
Kapelevich, Coey and Vielma, Coey, Kapelevich and Vielma (2022), and SOSTOOLS.

## License

MIT, see [LICENSE](LICENSE). SOSTOOLS is distributed separately under its own license (GPL) and is
not included here.
