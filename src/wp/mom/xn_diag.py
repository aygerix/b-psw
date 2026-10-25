#!/usr/bin/env python3
"""xn line: diagnostics at the least-squares minimiser of the order-N truncation (where is it, what is violated)."""
import sys
import numpy as np
from scipy.optimize import least_squares
from xn_core import E_both
from xn_feas import unpack, resid, margins

k, lin, N, ns, seed = map(int, sys.argv[1:6])
m = k - lin
rng = np.random.default_rng(seed)
best = None
for it in range(ns):
    p0 = np.concatenate([rng.normal(0, 1.5, m), rng.normal(0, 1.0, m)])
    with np.errstate(all="ignore"):
        sol = least_squares(resid, p0, args=(k, lin, N), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
    if np.all(np.isfinite(sol.fun)) and (best is None or sol.cost < best.cost):
        best = sol
u, v = unpack(best.x)
o = np.argsort(v)
print(f"k={k} lin={lin} N={N} best cost {best.cost:.6e}")
print("residual vector (scaled by 2^(k-1)):", np.array2string(best.fun, precision=6, max_line_width=250))
print("u =", np.array2string(u[o], precision=8, max_line_width=250))
print("v =", np.array2string(v[o], precision=8, max_line_width=250))
print("margins sqrt(8v)-(|u-1/2|+1/2):", np.array2string((np.sqrt(8 * v) - (np.abs(u - 0.5) + 0.5))[o], precision=3, max_line_width=250))
