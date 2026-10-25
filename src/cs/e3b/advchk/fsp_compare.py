#!/usr/bin/env python3
import os
"""src/cs/e3b/advchk/fsp_compare.py -- compare FSP point sets: E3b files (t, r_2..r_5) vs fsp_tu.py files (T, u_2..u_5).
Permutation-invariant family coordinates: (T, e_1(u), e_2(u), e_3(u), e_4(u)); clustering with relative tolerance tol.
Also re-polishes every point of the fsp_tu set with Newton (its own model) and reports the post-polish spread.
Usage: fsp_compare.py TYPE mine.txt theirs1.txt [theirs2.txt ...]"""
import sys, numpy as np
sys.path.insert(0, os.path.expanduser('~/b-psw/src/cs/e3b/advchk'))
import fsp_tu

def inv(T, u):
    e = np.poly(u)            # monic poly coefficients: [1, -e1, e2, -e3, e4]
    return np.array([T, -e[1], e[2], -e[3], e[4]])

def cluster(vs, tol):
    reps = []
    for v in vs:
        if not any(np.linalg.norm(v - w) <= tol * max(1, np.linalg.norm(w)) for w in reps):
            reps.append(v)
    return reps

typ = sys.argv[1]
mine_rows = [eval(l) for l in open(sys.argv[2])]
mine = []
for gd, z, near, sv in mine_rows:
    z = np.array(z)
    Z, nf = fsp_tu.newton(z[None].copy(), typ, iters=30)
    mine.append((inv(Z[0, 0], Z[0, 1:5]), nf[0], np.linalg.norm(Z[0] - z), sv))
theirs = []; drop = {'t~0': 0, 'nonconv': 0, 'deg': 0, 'dim>0': 0}
for fn in sys.argv[3:]:
    for l in open(fn):
        row = eval(l); zs = np.array(row[2])
        if row[0] != 0:
            drop['dim>0'] += 1
        t = zs[0]; r = zs[1:5]
        if abs(t) < 1e-3:
            drop['t~0'] += 1; continue
        z = np.concatenate([[t**2], (t - r) / t] + ([[zs[5] / t]] if typ != '1' else []))
        Z, nf = fsp_tu.newton(z[None].copy(), typ, iters=30)
        if not nf[0] < 1e-11:
            drop['nonconv'] += 1; continue
        if not fsp_tu.nondeg(Z[0]):
            drop['deg'] += 1; continue
        theirs.append(inv(Z[0, 0], Z[0, 1:5]))
print("their rows dropped:", drop)
print(f"mine: {len(mine)} rows; max polish move {max(m[2] for m in mine):.2e}; max residual {max(m[1] for m in mine):.2e}")
for tol in (1e-8, 1e-6, 1e-4, 1e-3):
    cm = cluster([m[0] for m in mine], tol); ct = cluster(theirs, tol)
    both = sum(1 for w in ct if any(np.linalg.norm(v - w) <= tol * max(1, np.linalg.norm(w)) for v in cm))
    print(f"tol {tol:g}: mine distinct {len(cm)}, theirs distinct {len(ct)} (rows {len(theirs)}), common {both}, "
          f"Lincoln-Petersen {len(cm) * len(ct) / max(both, 1):.1f}")
# conjugation-closure check of my set: is the conjugate of each of my points also in my set?
cm = cluster([m[0] for m in mine], 1e-6)
miss = sum(1 for v in cm if not any(np.linalg.norm(np.conj(v) - w) <= 1e-6 * max(1, np.linalg.norm(w)) for w in cm))
print(f"my points whose complex conjugate is NOT in my set: {miss} of {len(cm)}")
