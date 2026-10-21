#!/usr/bin/env python3
"""Apply the three Octave compatibility edits of README_OCTAVE.md to a copy of SOSTOOLS 4.

    python3 octave/patch_sostools.py /path/to/SOSTOOLS

Only needed to run SOSTOOLS under GNU Octave; MATLAB users need none of this. The edits change
no behaviour under MATLAB. The script is idempotent and reports each edit.
"""
import os
import sys

EDITS = [
    ("dpvar/@dpvar/compress.m",
     "C_new = C(Crows_new,Ccols_new);",
     "C_new = C(full(logical(Crows_new)),full(logical(Ccols_new)));"),
    ("dpvar/@dpvar/compress.m",
     "degmat_new = degmat(degrows,:);",
     "degmat_new = degmat(full(logical(degrows)),:);"),
    ("dpvar/@dpvar/private/DPVuniqueterm.m",
     "reshape(rvec + (0:nut:ncols_new-nut),ncols,1)",
     "reshape(full(rvec) + (0:nut:ncols_new-nut),ncols,1)"),
    ("dpvar/@dpvar/private/DPVuniqueterm.m",
     "reshape(sortidx + (0:nt:ncols-nt),ncols,1)",
     "reshape(full(sortidx) + (0:nt:ncols-nt),ncols,1)"),
    ("dpvar/@dpvar/diff.m",
     "y.*Jcoef",
     "y*spdiags(Jcoef(:),0,numel(Jcoef),numel(Jcoef))"),
]


def main(root):
    status = 0
    for rel, old, new in EDITS:
        path = os.path.join(root, rel)
        with open(path, newline="") as f:
            text = f.read()
        if new in text:
            print(f"already patched : {rel}")
            continue
        if old not in text:
            print(f"NOT FOUND       : {rel}: {old}")
            status = 1
            continue
        with open(path, "w", newline="") as f:
            f.write(text.replace(old, new))
        print(f"patched         : {rel}")
    return status


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1]))
