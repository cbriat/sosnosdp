# Running SOSNOSDP under GNU Octave

The package itself is plain MATLAB code that also runs in Octave (tested with Octave 8.4). SOSTOOLS 4,
however, relies on three MATLAB behaviours that Octave lacks. To use SOSTOOLS in Octave, edit a local
copy of SOSTOOLS as follows (MATLAB users need none of this):

1. `dpvar/@dpvar/compress.m`: Octave treats sparse logical indices as numeric. Replace
   `C_new = C(Crows_new,Ccols_new);` by `C_new = C(full(logical(Crows_new)),full(logical(Ccols_new)));`
   and `degmat_new = degmat(degrows,:);` by `degmat_new = degmat(full(logical(degrows)),:);`.
   Without this, every stored coefficient is silently wrong.
2. `dpvar/@dpvar/private/DPVuniqueterm.m`: Octave does not broadcast sparse arrays. Wrap `rvec` and
   `sortidx` in `full(...)` in the two `reshape(... + (0:...),ncols,1)` lines.
3. `dpvar/@dpvar/diff.m`: replace `y.*Jcoef` by `y*spdiags(Jcoef(:),0,numel(Jcoef),numel(Jcoef))`.

The script `octave/patch_sostools.py` applies these edits (`python3 octave/patch_sostools.py /path/to/SOSTOOLS`).

Also put `octave/contains.m` (a minimal replacement for MATLAB's `contains`) on the path.

Known Octave limitations of SOSTOOLS itself: `polynomial * dpvar` products fail (write `dpvar * polynomial`,
e.g. `s * (x*(2 - x))`), and `sosmatrixineq(prog, M)` in its default form fails, so the test suite uses an
equivalent reconstruction (`tests/quadratic_matrix_ineq.m`). SOSTOOLS also needs an SDP solver on the
path to build programs; `setup_nosdp` adds a placeholder when none is installed.

Run the tests with

    SOSTOOLS_DIR = '/path/to/patched/SOSTOOLS'; run('octave/run_tests_octave.m')
