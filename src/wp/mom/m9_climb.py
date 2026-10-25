#!/usr/bin/env python3
"""m9 line: climb the Taylor hierarchy (levels J = J0, J0+1, ...; orders 0..J at t = -1 for both halves) from a
given solution, by least squares continuation with random perturbations at each level.  Float64; a level counts as
reached if max residual / c0 < tol, min domain margin > 1e-6 and min member distance > 1e-4.
Usage: m9_climb.py k lin J0 Jmax tries seed 'u,..' 'v,..' [tol]
"""
import sys
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, '.')
from m9_hier import resid, unpack, taylor_resid
k, lin, J0, Jmax, tries, seed = map(int, sys.argv[1:7])
U = np.array([float(x) for x in sys.argv[7].split(',')]); V = np.array([float(x) for x in sys.argv[8].split(',')])
tol = float(sys.argv[9]) if len(sys.argv) > 9 else 1e-7
m = k - lin
rng = np.random.default_rng(seed)
r0 = np.sqrt(8 * V) - .5
q = np.concatenate([np.log(r0), np.arctanh(np.clip((U - .5) / r0, -.999999, .999999))])
np.set_printoptions(precision=16, linewidth=300)
for J in range(J0, Jmax + 1):
    got = None
    for t in range(tries):
        q0 = q + (0 if t == 0 else rng.normal(0, 0.05 * t, 2 * m))
        try:
            s = least_squares(resid, q0, args=(k, lin, J), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=8000)
        except Exception:
            continue
        u, v = unpack(s.x, m)
        f = np.abs(s.fun).max()
        mg = (np.sqrt(8 * v) - .5 - np.abs(u - .5)).min()
        dm = min(np.hypot(u[i] - u[j], v[i] - v[j]) for i in range(m) for j in range(i + 1, m))
        print(f"level {J} try {t}: maxres {f:.2e} margin {mg:.3e} dmin {dm:.3e} S {v.sum():.3f} vmax {v.max():.3f}", flush=True)
        if f < tol and mg > 1e-6 and dm > 1e-4:
            got = s.x
            break
    if got is None:
        print(f"LEVEL {J}: NOT reached (best attempts above)", flush=True)
        break
    q = got
    u, v = unpack(q, m)
    print(f"LEVEL {J}: reached", flush=True)
    print("   u =", repr(u)); print("   v =", repr(v), flush=True)
