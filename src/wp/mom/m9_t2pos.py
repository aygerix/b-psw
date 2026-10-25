#!/usr/bin/env python3
"""m9 line: search the level-1 set (four boundary equations) for points with 2^k T2+ = theta > 0, including
LARGE members (initial rho ~ N(mu, 1.5)).  5 equations, least squares, exact domain parametrisation.
Usage: m9_t2pos.py k lin theta mu nstarts seed
"""
import sys
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, '.')
from m9_feas import resid as resid1
from m9_ord2 import quantities
k, lin = int(sys.argv[1]), int(sys.argv[2])
th, mu = float(sys.argv[3]), float(sys.argv[4])
ns, seed = int(sys.argv[5]), int(sys.argv[6])
m = k - lin
rng = np.random.default_rng(seed)
def r5(q):
    u, v, a, b, es, ea = quantities(q, k, lin)
    return np.append(resid1(q, k, lin), ((a + b) / 2 - th) / 2.0 ** (k - 1))
nsol = 0
for it in range(ns):
    q0 = np.concatenate([rng.normal(mu, 1.5, m), rng.normal(0, 1.0, m)])
    try:
        s = least_squares(r5, q0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=5000)
    except Exception:
        continue
    u, v, a, b, es, ea = quantities(s.x, k, lin)
    f = np.abs(s.fun).max()
    mg = (np.sqrt(8 * v) - .5 - np.abs(u - .5)).min()
    dm = min(np.hypot(u[i] - u[j], v[i] - v[j]) for i in range(m) for j in range(i + 1, m))
    ok = f < 1e-7 and mg > 1e-6 and dm > 1e-4
    nsol += ok
    print(f"start {it}: res {f:.2e} S {v.sum():.3f} vmax {v.max():.3f} 2^kT2+ {(a + b) / 2:.4e} margin {mg:.2e} dmin {dm:.2e} {'SOL' if ok else 'fail'}", flush=True)
    if ok:
        print("   u =", repr(u)); print("   v =", repr(v))
print("solutions with 2^k T2+ =", th, ":", nsol, "of", ns)
