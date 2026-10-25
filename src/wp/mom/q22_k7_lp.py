"""For fixed (r1, r2), look for weights alpha in R^8 with sum_i alpha_i B_i(w) >= 1 for all sampled w > 0 (7 members).
If such alpha exists on the whole domain, the 8 boundary equations have no solution with c_j > 0 at that shape.
Usage: python3 q22_k7_lp.py ngrid nsamples"""
import sys, numpy as np
from scipy.optimize import linprog
exec(open('q22_k7_feas.py').read().split('starts, seed')[0])
ng, ns = int(sys.argv[1]), int(sys.argv[2]); rng = np.random.default_rng(11)
grid = np.linspace(-0.95, 0.95, ng)
bad = []; ok = 0
for r1 in grid:
    for r2 in grid:
        if r2 <= r1 + 1e-9: continue
        W = np.exp(rng.uniform(-7, 5, (ns, k)))
        Bm = np.array([[eval(c, {}, env(r1, r2, w)) for c in code] for w in W])
        Bm = Bm / np.maximum(1.0, np.abs(Bm).max(axis=1, keepdims=True))      # scale rows
        res = linprog(np.zeros(8), A_ub=-Bm, b_ub=-1e-6 * np.ones(ns), bounds=[(-1e6, 1e6)] * 8, method='highs')
        if res.status == 0: ok += 1
        else: bad.append((round(r1, 3), round(r2, 3)))
print(f"shapes tested {ok + len(bad)}: separating weights found {ok}; not found {len(bad)}")
print("not found at:", bad[:40])
