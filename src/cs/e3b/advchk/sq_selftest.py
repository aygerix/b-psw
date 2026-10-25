#!/usr/bin/env python3
import os
"""src/cs/e3b/advchk/sq_selftest.py -- does sq_enum.Sys represent the known symmetric (2,2,2,0,0) points of
out/sym_22200.txt (structure 1|1|2|2|3|3)?  Residual of the full system at the known point and basin of the SQUARED
system from perturbations.  Usage: sq_selftest.py"""
import sys, numpy as np
sys.path.insert(0, os.path.expanduser('~/b-psw/src/cs/e3b/advchk'))
sys.path.insert(0, os.path.expanduser('~/b-psw/src/cs/e3b'))
import sq_enum
from refine import lucas_table
rng = np.random.default_rng(7)
S = sq_enum.Sys("1|1|2|2|3|3", rng=rng)
rows = [eval(l) for l in open(os.path.expanduser('~/b-psw/src/cs/e3b/out/sym_22200.txt'))][:6]
for d, gd, clr, pts in rows:
    clr = np.array(clr)
    tab = lucas_table(clr, 5)
    lp = [[t for (t, r) in row if r < 1e-6] for row in tab]
    z = np.concatenate([clr, [lp[0][0], -lp[0][0], lp[1][0], lp[1][1], lp[2][0], lp[2][1]]])
    r0 = np.linalg.norm(S.full(z[None])[0])
    line = f"full residual {r0:.1e}; Lucas recount {[len(p) for p in lp]};"
    for eps in (1e-3, 1e-2, 1e-1, 0.3):
        Z = z[None] + eps * (rng.normal(size=(100, S.n)) + 1j * rng.normal(size=(100, S.n)))
        with np.errstate(all='ignore'):
            Z, nf = sq_enum.newton(S, Z.copy())
        back = np.sum((nf < 1e-11) & (np.linalg.norm(Z - z[None], axis=1) < 1e-5))
        line += f" eps={eps:g}: {back}/100"
    print(line, flush=True)
