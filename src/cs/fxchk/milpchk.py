"""Independent exact solver for M(lambda) = max_C [N(C) - lambda bits(C)] on a bbc instance file (adversarial FX check).
MILP (HiGHS via scipy.optimize.milp): x_p, y_item continuous in [0,1]; one binary z_q per row needed on BOTH sides by some
prime (F side allowed iff z_q = 1).  With z fixed the remaining system is a closure polytope (difference constraints,
totally unimodular), so only z needs integrality.  The incumbent is re-evaluated exactly (own code) from the primes.
Usage: python3 milpchk.py INST lam1[,lam2] [time_limit]"""
import sys, os, time, math
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from instlib import read_inst, item_q
fn = sys.argv[1]; lams = [float(x) for x in sys.argv[2].split(',')]; TL = float(sys.argv[3]) if len(sys.argv) > 3 else 3000
nP, nI, nR, items, prim = read_inst(fn); qs = item_q(items)
cost = [math.log2(qs[j] - 1) if items[j][2] == 1 else math.log2(qs[j]) for j in range(nI)]
# rows needed on both sides
sides = {}
for p, its in prim:
    for j in its: sides.setdefault(items[j][0], set()).add(items[j][1])
conf = sorted(r for r, s in sides.items() if len(s) == 2)
zid = {r: nP + nI + k for k, r in enumerate(conf)}
nv = nP + nI + len(conf)
rows, cols, vals, ub = [], [], [], []
k = 0
for i, (p, its) in enumerate(prim):
    for j in its: rows += [k, k]; cols += [i, nP + j]; vals += [1, -1]; ub.append(0); k += 1
for j, (r, s, e, par, c) in enumerate(items):
    if par >= 0: rows += [k, k]; cols += [nP + j, nP + par]; vals += [1, -1]; ub.append(0); k += 1
    if e == 1 and r in zid:
        if s == 0: rows += [k, k]; cols += [nP + j, zid[r]]; vals += [1, -1]; ub.append(0); k += 1
        else: rows += [k, k]; cols += [nP + j, zid[r]]; vals += [1, 1]; ub.append(1); k += 1
A = coo_matrix((vals, (rows, cols)), shape=(k, nv)).tocsr()
integ = np.zeros(nv); integ[nP + nI:] = 1
def exact(sel_primes):
    used = set()
    for i in sel_primes:
        for j in prim[i][1]:
            while j >= 0 and j not in used: used.add(j); j = items[j][3]
    sd = {}
    ok = True
    for j in used:
        if items[j][2] == 1:
            if items[j][0] in sd and sd[items[j][0]] != items[j][1]: ok = False
            sd[items[j][0]] = items[j][1]
    return len(sel_primes), sum(cost[j] for j in used), ok
print(f"{os.path.basename(fn)}: primes {nP} items {nI} rows {nR} conflict rows {len(conf)}")
for lam in lams:
    c = np.zeros(nv); c[:nP] = -1; c[nP:nP + nI] = lam * np.array(cost)
    t0 = time.time()
    res = milp(c, constraints=[LinearConstraint(A, -np.inf, np.array(ub, dtype=float))], integrality=integ, bounds=Bounds(0, 1),
               options=dict(time_limit=TL, mip_rel_gap=0.0, disp=False))
    dt = time.time() - t0
    if res.x is None: print(f"lam {lam}: status {res.status} {res.message} time {dt:.1f}"); continue
    sel = [i for i in range(nP) if res.x[i] > 0.5]
    N, bits, ok = exact(sel)
    db = getattr(res, 'mip_dual_bound', None)
    print(f"lam {lam:.4f}: status {res.status} opt {-res.fun:.6f} dual_bound {(-db if db is not None else float('nan')):.6f} "
          f"incumbent exact N {N} bits {bits:.4f} value {N - lam * bits:.6f} sides_ok {ok} time {dt:.1f} nodes {getattr(res, 'mip_node_count', '?')}", flush=True)
