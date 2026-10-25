#!/usr/bin/env python3
"""xn line: a priori bound on S = sum 1/a_i from the expansion at t = -1 (k = 9).

For a family V = Vpoly, so V's Taylor coefficients at t = -1, V_n = [tau^n] (1 - tau)^8 G(tau)  (G-side: polynomials
in (a_i, c_i), bounded on the domain), determine all coefficients of V.  With Vpoly = t^8 + 0 t^7 + 0 t^6 + S t^5 + ...:
   1 = V_8,  0 = V_7 + 8 V_8,  0 = V_6 + 7 V_7 + 28 V_8,  S = V_5 + 6 V_6 + 21 V_7 + 56 V_8.
So S = SG := V_5 + 6 V_6 + 21 V_7 + 56 V_8 (G-side, no family assumption needed to evaluate it).
This script maximises SG numerically over the real domain (COMPUTED); a rigorous bound is separate.
Usage: xn_Sbound.py lin nstarts seed
"""
import sys
import numpy as np
from scipy.optimize import minimize
from xn_core import G_tau

lin, ns, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
k = 9
m = k - lin
rng = np.random.default_rng(seed)


def unpack(p):
    # a in (0, 32); c with c^2 < 8a and (a - c)^2 < 8a  <=>  c in (a - sqrt(8a), sqrt(8a)) intersect (-sqrt(8a), ...)
    a = 32 / (1 + np.exp(-p[:m]))
    lo = np.maximum(-np.sqrt(8 * a), a - np.sqrt(8 * a))
    hi = np.minimum(np.sqrt(8 * a), a + np.sqrt(8 * a))
    c = lo + (hi - lo) / (1 + np.exp(-p[m:]))
    return a, c


def SG(p, sign=1):
    a, c = unpack(p)
    Gt = G_tau(list(a), list(c), lin, k, 9)
    val = Gt[5] + 6 * Gt[6] + 21 * Gt[7] + 56 * Gt[8]
    return -sign * val


best = -np.inf
bestmin = np.inf
for it in range(ns):
    p0 = rng.normal(0, 2, 2 * m)
    r = minimize(SG, p0, method="Nelder-Mead", options={"maxiter": 20000, "xatol": 1e-9, "fatol": 1e-12})
    r = minimize(SG, r.x, method="BFGS")
    best = max(best, -r.fun)
    r2 = minimize(lambda p: SG(p, -1), p0, method="Nelder-Mead", options={"maxiter": 20000})
    bestmin = min(bestmin, r2.fun)
    a, c = unpack(r.x)
    print(f"start {it}: max SG {-r.fun:.6g}  (a = {np.array2string(np.sort(a), precision=3)})   min SG {r2.fun:.6g}", flush=True)
print(f"lin={lin}: sup SG found {best:.6g}, inf SG found {bestmin:.6g}")
