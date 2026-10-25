#!/usr/bin/env python3
"""m9 line: solve hierarchy level J (Taylor orders 0..J at t = -1, both halves) by least squares from a given start
(u, v) plus random perturbations; report solutions inside the domain.  Uses m9_hier.resid (stable series version).
Usage: m9_lev2.py k lin J ntrials seed 'u0,u1,...' 'v0,v1,...'
"""
import sys
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, '.')
from m9_hier import resid, unpack, taylor_resid

k, lin, J, nt, seed = map(int, sys.argv[1:6])
U = np.array([float(x) for x in sys.argv[6].split(',')])
V = np.array([float(x) for x in sys.argv[7].split(',')])
m = k - lin
rng = np.random.default_rng(seed)
r0 = np.sqrt(8 * V) - .5
q0 = np.concatenate([np.log(r0), np.arctanh(np.clip((U - .5) / r0, -.999999, .999999))])
for it in range(nt):
    q = q0 + (0 if it == 0 else rng.normal(0, 0.1, 2 * m))
    s = least_squares(resid, q, args=(k, lin, J), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=8000)
    u, v = unpack(s.x, m)
    f = np.abs(s.fun).max()
    mg = (np.sqrt(8 * v) - .5 - np.abs(u - .5)).min()
    dm = min(np.hypot(u[i] - u[j], v[i] - v[j]) for i in range(m) for j in range(i + 1, m))
    ok = f < 1e-7 and mg > 1e-6 and dm > 1e-4
    print(f"trial {it}: level {J} maxres {f:.2e} margin {mg:.3e} dmin {dm:.3e} {'SOL' if ok else 'fail'}", flush=True)
    if ok:
        np.set_printoptions(precision=16, linewidth=300)
        print("   u =", repr(u)); print("   v =", repr(v))
        a, _ = taylor_resid(u, v, k, lin, J + 3); b, _ = taylor_resid(1 - u, v, k, lin, J + 3)
        print("   next orders / c0:", a[J + 1:] / 2.0 ** (1 - k), b[J + 1:] / 2.0 ** (1 - k))
