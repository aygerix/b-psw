#!/usr/bin/env python3
"""src/cs/e3b/sym_mp.py -- E3b: high-precision refinement of REAL points from strat_lp output files (any structure,
built from the point's own Lucas points) and the rationality test of the family (mp_refine.invariants: c, l_j r_j,
r_i r_j all rational <=> defined over Q after x -> sqrt(d) x + h).  Distinct points only (up to x -> -x, relabelling).
Usage: sym_mp.py file [dps]"""
import sys, numpy as np, mpmath as mp
sys.path.insert(0, __file__.rsplit('/', 1)[0] if '/' in __file__ else '.')
from node_model import Node
from mp_refine import newton, invariants, mindeg
from refine import lucas_table
fn = sys.argv[1]; dps = int(sys.argv[2]) if len(sys.argv) > 2 else 60
seen = set()
for line in open(fn):
    d, gd, clr, pts = eval(line)
    clr = np.array(clr)
    if np.max(np.abs(clr.imag)) > 1e-8:
        continue
    key = tuple(np.round(sorted(np.abs(clr)), 4))
    if key in seen:
        continue
    seen.add(key)
    tab = lucas_table(clr, 5)
    lp = [[t for (t, r) in row if r < 1e-6] for row in tab]
    # node structure: one node per (member, Lucas point); member 1 uses (w, -w)
    nodes = []; pos = []
    for j in range(5):
        for t in lp[j]:
            nodes.append((j + 1,)); pos.append(t)
    model = Node(nodes)
    z = np.zeros(model.n, dtype=complex); z[:4] = clr[5:9]
    for i, t in enumerate(pos):
        z[4 + i] = t
    for j, ix in model.coef.items():
        z[ix] = clr[0] if j == 1 else clr[j - 1]
    zz, nf = newton(model, list(z), dps=dps)
    inv, c, lam, rr = invariants(model, zz)
    degs = {k: mindeg(v) for k, v in inv.items()}
    rat = all(isinstance(v, tuple) for v in degs.values())
    print(f"gcd {[len(p) for p in lp]} resid {mp.nstr(nf, 3)} rational_family={rat}  c = {mp.nstr(mp.re(c), 20)}")
    print("    degrees:", {k: (v if not isinstance(v, tuple) else 'Q') for k, v in degs.items()}, flush=True)
