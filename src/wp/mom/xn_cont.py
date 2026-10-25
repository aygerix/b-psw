#!/usr/bin/env python3
"""xn line: order-by-order continuation for the truncated Williams system (real cubic families).

Stage N = 2 (boundary equations): many random starts.  Every converged point is used as the start of stage
N + 1 (all 2(N+1) equations), and so on up to Nmax.  A stage "converges" if maxres < 1e-9 (scaled by 2^(k-1)).
At each stage also random restarts (nrand) are tried.  Reports the largest N reached and the points.
Usage: xn_cont.py k lin Nmax nstarts nrand seed
COMPUTED evidence only.
"""
import sys
import numpy as np
from scipy.optimize import least_squares
from xn_feas import unpack, resid, margins

k, lin, Nmax, ns, nrand, seed = map(int, sys.argv[1:7])
m = k - lin
rng = np.random.default_rng(seed)


def solve(p0, N, mx=4000):
    with np.errstate(all="ignore"):
        try:
            sol = least_squares(resid, p0, args=(k, lin, N), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=mx)
        except Exception:
            return None, np.inf
    if not np.all(np.isfinite(sol.fun)):
        return None, np.inf
    return sol.x, np.abs(sol.fun).max()


pts = []
for it in range(ns):
    p0 = np.concatenate([rng.normal(0, 1.5, m), rng.normal(0, 1.0, m)])
    x, f = solve(p0, 2)
    if f < 1e-9:
        u, v = unpack(x)
        mg, dm = margins(u, v)
        if dm > 1e-6 and mg > 1e-8:
            pts.append(x)
print(f"k={k} lin={lin}: N=2 solutions (interior, distinct): {len(pts)} of {ns}", flush=True)
for N in range(3, Nmax + 1):
    new, best = [], np.inf
    starts = pts + [np.concatenate([rng.normal(0, 1.5, m), rng.normal(0, 1.0, m)]) for _ in range(nrand)]
    for p0 in starts:
        x, f = solve(p0, N)
        best = min(best, f)
        if f < 1e-9:
            u, v = unpack(x)
            mg, dm = margins(u, v)
            if dm > 1e-6 and mg > 1e-8:
                new.append(x)
    print(f"N={N}: solutions {len(new)} from {len(starts)} starts; best maxres {best:.3e}", flush=True)
    for x in new[:3]:
        u, v = unpack(x)
        o = np.argsort(v)
        mg, dm = margins(u, v)
        print("   margin %.3e dmin %.3e" % (mg, dm))
        print("   u =", np.array2string(u[o], precision=10, max_line_width=300))
        print("   v =", np.array2string(v[o], precision=10, max_line_width=300))
    if not new:
        print(f"no solution at N={N}")
        break
    pts = new
