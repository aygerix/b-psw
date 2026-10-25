#!/usr/bin/env python3
"""src/cs/e3b/classify.py -- E3b: classify stratum solutions: x -> -x symmetric (r's and l's in mirror pairs, c any),
FSP (a Lucas point common to all 5 members), real/complex, and cluster distinct points (up to x -> -x and relabelling).
Usage: classify.py file k"""
import sys, numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0] if '/' in __file__ else '.')
from refine import lucas_table
from collections import Counter
fn, k = sys.argv[1], int(sys.argv[2])
rows = [eval(l) for l in open(fn)]
cl = Counter(); reps = {}
for d, gd, clr, pts in rows:
    clr = np.array(clr); c = clr[0]; lam = clr[1:k]; rr = clr[k:]
    tab = lucas_table(clr, k)
    lp = [[t for (t, r) in row if r < 1e-6] for row in tab]
    gd2 = [len(p) for p in lp]
    # symmetric: there is a fixed-point-free pairing of cubic members with l' = -l, r' = -r
    sym = False
    idx = list(range(k - 1)); used = set(); ok = True
    for i in idx:
        if i in used: continue
        mate = [j for j in idx if j != i and j not in used and abs(lam[j] + lam[i]) < 1e-6 * max(1, abs(lam[i])) and abs(rr[j] + rr[i]) < 1e-6 * max(1, abs(rr[i]))]
        if not mate: ok = False; break
        used |= {i, mate[0]}
    sym = ok
    allpts = [t for p in lp for t in p]
    fsp = any(all(any(abs(t - u) < 1e-6 for u in p) for p in lp) for t in lp[0]) if lp[0] else False
    real = np.max(np.abs(clr.imag)) < 1e-8
    inv = tuple(np.round(sorted([abs(c)] + list(np.abs(rr)) + list(np.abs(lam))), 5))
    key = (tuple(sorted(gd2[1:], reverse=True)), gd2[0], 'sym' if sym else 'nosym', 'fsp' if fsp else '-', 'real' if real else 'cplx')
    cl[(key, inv)] += 1
    reps.setdefault((key, inv), clr)
groups = Counter()
for (key, inv), v in cl.items():
    groups[key] += 1
print("distinct points (up to x->-x, relabel) by class:")
for key, v in sorted(groups.items()):
    print("  ", key, ":", v, "distinct points")
for (key, inv), clr in reps.items():
    if key[4] == 'real':
        print("  real rep", key, np.round(clr.real, 6).tolist())
