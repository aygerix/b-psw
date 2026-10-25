"""Independent exact solve of M(lambda) = max_C [N(C) - lambda bits(C)] for one (D, 2-adic class) by MILP (scipy milp /
HiGHS), built directly from the record files (no instance file, no min cut, no deletion lemma, no lagr_bb code).
Variables: y_{q,s,e} binary (row q on side s with cap >= e), x_p in [0,1].  Constraints: x_p <= y_{q,s,a} for every
requirement q^a || o_s(p); y_{q,s,e} <= y_{q,s,e-1}; y_{q,0,1} + y_{q,1,1} <= 1.  Objective max sum x - lambda sum c y,
c(q,1) = log2(q-1), c(q,e>1) = log2 q.  Only primes whose odd orders are T0-smooth are kept (rows <= T0).
Usage: python3 aud_milp.py D j1 j2 T0 lam1[,lam2] timelimit recfile1 [recfile2 ...]"""
import sys, math, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix
D, j1, j2, T0 = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
lams = [float(x) for x in sys.argv[5].split(',')]; TL = float(sys.argv[6]); files = sys.argv[7:]
PR = [q for q in range(3, 201) if all(q % d for d in range(2, int(q ** .5) + 1))]
def fac(x):
    f = {}
    for q in PR:
        while x % q == 0: x //= q; f[q] = f.get(q, 0) + 1
    return f if x == 1 else None
pool = {}
for fn in files:
    for l in open(fn):
        t = tuple(map(int, l.split()))
        if len(t) != 7 or t[0] != D or (t[3], t[4]) != (j1, j2): continue
        pool[t[1]] = (t[5], t[6])
reqs = {}
for p, (o1, o2) in pool.items():
    f1, f2 = fac(o1), fac(o2)
    if any(q > T0 for q in f1) or any(q > T0 for q in f2): continue
    reqs[p] = [(q, 0, a) for q, a in f1.items()] + [(q, 1, a) for q, a in f2.items()]
P = sorted(reqs); nP = len(P)
emax = {}
for p in P:
    for q, s, a in reqs[p]: emax[(q, s)] = max(emax.get((q, s), 0), a)
items = [(q, s, e) for (q, s), E in sorted(emax.items()) for e in range(1, E + 1)]
iid = {it: nP + k for k, it in enumerate(items)}
nv = nP + len(items)
cost = np.array([math.log2(q - 1) if e == 1 else math.log2(q) for q, s, e in items])
rr, cc, vv, ub = [], [], [], []; m = 0
for k, p in enumerate(P):
    for q, s, a in reqs[p]:
        rr += [m, m]; cc += [k, iid[(q, s, a)]]; vv += [1, -1]; ub.append(0); m += 1
for (q, s, e) in items:
    if e > 1: rr += [m, m]; cc += [iid[(q, s, e)], iid[(q, s, e - 1)]]; vv += [1, -1]; ub.append(0); m += 1
for q in sorted({q for q, s in emax}):
    if (q, 0) in emax and (q, 1) in emax:
        rr += [m, m]; cc += [iid[(q, 0, 1)], iid[(q, 1, 1)]]; vv += [1, 1]; ub.append(1); m += 1
A = coo_matrix((vv, (rr, cc)), shape=(m, nv)).tocsr()
integ = np.concatenate([np.zeros(nP), np.ones(len(items))])
print(f"D = {D}, class ({j1},{j2}), T0 = {T0}: primes {nP}, items {len(items)}, rows {len({q for q, s in emax})}, "
      f"constraints {m}", flush=True)
for lam in lams:
    c = np.concatenate([-np.ones(nP), lam * cost])
    t0 = time.time()
    r = milp(c, constraints=LinearConstraint(A, -np.inf, np.array(ub)), integrality=integ, bounds=Bounds(0, 1),
             options=dict(time_limit=TL, disp=False, mip_rel_gap=0))
    val = None if r.x is None else -r.fun
    dual = getattr(r, 'mip_dual_bound', None)
    dual = None if dual is None else -dual
    N = bits = None
    if r.x is not None:
        y = r.x[nP:] > 0.5; N = int(round(r.x[:nP].sum())); bits = float(cost[y].sum())
        sel = [P[k] for k in range(nP) if r.x[k] > 0.5]
    print(f"lam {lam:.4f}: status {r.status} ({r.message.split('.')[0]}), primal M = {val}, dual bound = {dual}, "
          f"at (bits {bits}, N {N}), {time.time() - t0:.1f} s; primes {sel[:20] if r.x is not None else None}", flush=True)
