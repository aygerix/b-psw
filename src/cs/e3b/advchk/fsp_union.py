#!/usr/bin/env python3
import os
"""src/cs/e3b/advchk/fsp_union.py -- population estimate for an FSP point type from several fsp_tu.py runs plus E3b runs.
Family invariant: (T, e_1..e_4 of u).  Prints per-run distinct counts, pairwise Lincoln-Petersen between my runs, the
union, and the union closed under complex conjugation (the true set is conjugation-closed: real equations).
Usage: fsp_union.py TYPE mine1.txt,mine2.txt[,...] [theirs1.txt,theirs2.txt,...]"""
import sys, itertools, numpy as np
sys.path.insert(0, os.path.expanduser('~/b-psw/src/cs/e3b/advchk'))
import fsp_tu

def inv(T, u):
    e = np.poly(u)
    return np.array([T, -e[1], e[2], -e[3], e[4]])

def same(v, w, tol=1e-7):
    return np.linalg.norm(v - w) <= tol * max(1, np.linalg.norm(w))

def dedup(vs):
    out = []
    for v in vs:
        if not any(same(v, w) for w in out):
            out.append(v)
    return out

typ = sys.argv[1]
runs = []
for fn in sys.argv[2].split(','):
    runs.append(dedup([inv(np.array(eval(l)[1])[0], np.array(eval(l)[1])[1:5]) for l in open(fn)]))
theirs = []
if len(sys.argv) > 3:
    for fn in sys.argv[3].split(','):
        for l in open(fn):
            row = eval(l); zs = np.array(row[2]); t = zs[0]; r = zs[1:5]
            if abs(t) < 1e-3:
                continue
            z = np.concatenate([[t**2], (t - r) / t] + ([[zs[5] / t]] if typ != '1' else []))
            Z, nf = fsp_tu.newton(z[None].copy(), typ, iters=30)
            if nf[0] < 1e-11 and fsp_tu.nondeg(Z[0]):
                theirs.append(inv(Z[0, 0], Z[0, 1:5]))
    theirs = dedup(theirs)
for i, r in enumerate(runs):
    print(f"my run {i}: {len(r)} distinct")
print(f"E3b union: {len(theirs)} distinct")
for i, j in itertools.combinations(range(len(runs)), 2):
    both = sum(1 for v in runs[i] if any(same(v, w) for w in runs[j]))
    print(f"runs {i},{j}: common {both}, Lincoln-Petersen {len(runs[i]) * len(runs[j]) / max(both, 1):.0f}")
allmine = dedup([v for r in runs for v in r])
if theirs:
    both = sum(1 for v in theirs if any(same(v, w) for w in allmine))
    print(f"my union {len(allmine)} vs E3b union {len(theirs)}: common {both}, Lincoln-Petersen {len(allmine) * len(theirs) / max(both, 1):.0f}")
U = dedup(allmine + theirs)
Uc = dedup(U + [np.conj(v) for v in U])
nreal = sum(1 for v in Uc if np.max(np.abs(v.imag)) < 1e-8)
print(f"union {len(U)}; union closed under conjugation {len(Uc)} (real {nreal})")
