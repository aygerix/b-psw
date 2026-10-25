#!/usr/bin/env python3
"""src/cs/e3b/fsp_dimchk.py -- E3b: slicing test of the local dimension of FSP-model solutions.
For each distinct point in a fsp_num output file: add d random affine slices through a perturbed copy of the point and
re-converge (LM) from the point; the component has dim >= d iff this succeeds (3 trials).  Usage: fsp_dimchk.py extra file"""
import sys, numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0] if '/' in __file__ else '.')
from fsp_num import FSP, to_clr
from strat_lp import lm_batch
from refine import lucas_table
from collections import Counter
extra = [e for e in sys.argv[1].split(',') if e]
rows = [eval(l) for l in open(sys.argv[2])]
rng = np.random.default_rng(3)
res = Counter()
for d0, gd, zs in rows:
    zs = np.array(zs)
    best = 0
    for d in (1, 2):
        ok = 0
        for trial in range(3):
            m = FSP(extra, nslice=d, rng=rng)
            zp = zs + 1e-3 * (rng.normal(size=zs.size) + 1j * rng.normal(size=zs.size))
            m.s0 = m.S @ zp
            zz, nf = lm_batch(m, zs[None].copy(), iters=300)
            if nf[0] < 1e-11 and np.linalg.norm(zz[0] - zs) < 0.05:
                ok += 1
        if ok == 0:
            break
        best = d
    tab = lucas_table(to_clr(zs), 5)
    gd2 = [sum(1 for (t, r) in row if r < 1e-8) for row in tab]
    res[(best, tuple(gd2))] += 1
    if best >= 1:
        print("POSITIVE-DIM candidate:", best, gd2, np.round(zs, 6).tolist())
for k, v in sorted(res.items()):
    print("  slice-dim", k[0], "gcd", list(k[1]), ":", v)
