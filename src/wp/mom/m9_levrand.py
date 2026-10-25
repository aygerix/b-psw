#!/usr/bin/env python3
"""m9 line: hierarchy level J (orders 0..J at t = -1, both halves) from RANDOM starts with large members
(rho ~ N(mu, 1.5)); counts solutions inside the domain.  Usage: m9_levrand.py k lin J mu nstarts seed [tol]"""
import sys
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, '.')
from m9_hier import resid, unpack
k, lin, J = map(int, sys.argv[1:4]); mu = float(sys.argv[4]); ns, seed = map(int, sys.argv[5:7])
tol = float(sys.argv[7]) if len(sys.argv) > 7 else 1e-6
m = k - lin
rng = np.random.default_rng(seed)
np.set_printoptions(precision=16, linewidth=300)
nok = 0
for it in range(ns):
    q0 = np.concatenate([rng.normal(mu, 1.5, m), rng.normal(0, 1.0, m)])
    try:
        s = least_squares(resid, q0, args=(k, lin, J), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=6000)
    except Exception:
        continue
    u, v = unpack(s.x, m)
    f = np.abs(s.fun).max()
    mg = (np.sqrt(8 * v) - .5 - np.abs(u - .5)).min()
    dm = min(np.hypot(u[i] - u[j], v[i] - v[j]) for i in range(m) for j in range(i + 1, m))
    ok = f < tol and mg > 1e-6 and dm > 1e-4
    nok += ok
    print(f"start {it}: level {J} maxres {f:.2e} margin {mg:.2e} dmin {dm:.2e} S {v.sum():.2f} vmax {v.max():.2f} {'SOL' if ok else 'fail'}", flush=True)
    if ok:
        print("   u =", repr(u)); print("   v =", repr(v), flush=True)
print(f"level {J}: {nok} solutions of {ns}")
