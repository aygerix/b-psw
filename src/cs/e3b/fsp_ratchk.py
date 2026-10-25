#!/usr/bin/env python3
"""src/cs/e3b/fsp_ratchk.py -- E3b: for FSP-model points (t, r2..r5) list the invariants T = t^2, u_i = (t - r_i)/t
(a family over Q <=> all of them rational) for the real points, with best rational approximations (denominator <= 1e4).
Usage: fsp_ratchk.py file"""
import sys, numpy as np
from fractions import Fraction
rows = [eval(l) for l in open(sys.argv[1])]
nreal = 0
for d, gd, zs in rows:
    zs = np.array(zs); t = zs[0]; rr = zs[1:5]
    T = t**2; u = (t - rr) / t
    inv = np.array([T] + list(u))
    if np.max(np.abs(inv.imag)) > 1e-8:
        continue
    nreal += 1
    inv = inv.real
    approx = [Fraction(v).limit_denominator(10**4) for v in inv]
    err = max(abs(float(a) - v) for a, v in zip(approx, inv))
    print(gd, "T,u =", np.round(inv, 8).tolist(), " best rationals", [str(a) for a in approx], f"err {err:.1e}")
print("real points:", nreal, "of", len(rows))
