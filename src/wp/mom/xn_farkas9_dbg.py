#!/usr/bin/env python3
"""debug helper: evaluate the rounded certificate Psi on a fine (r, theta) grid and report the minimum"""
import sys
import numpy as np
from fractions import Fraction as Fr
import xn_farkas9 as X

lin, Sa, Sb = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
flags = sys.argv[4].split("+") if sys.argv[4] != "none" else []
nv = int(sys.argv[5]) if len(sys.argv) > 5 else 120
nw = int(sys.argv[6]) if len(sys.argv) > 6 else 9
m = 9 - lin
Bh_list = None
t0_list = [Fr(-4) + Fr(7 * i, 35) for i in range(36)] if "POS" in flags else None
if "N3L" in flags:
    lo = 2 * m + Fr(m * m) / (2 * Fr(Sb).limit_denominator(10**6)) + 2 * lin
    Bh_list = [lo + (18 * m + 2 * lin - lo) * Fr(i, 39) for i in range(40)]
m, eqs, ineqs, names, res = X.farkas_lp(lin, Sa, Sb, flags, nv, nw, Bh_list, t0_list)
print("status", res.status, "eps", res.x[-1] if res.status == 0 else None)
ne, ni = len(eqs), len(ineqs)
y, z = res.x[:ne], res.x[ne:ne + ni]
print("y =", y); print("z nonzero:", [(names[j], z[j]) for j in range(ni) if z[j] > 1e-9][:10])
for S in (Sa, Sb):
    vmax = Sb - (m - 1) / 32
    R = np.sqrt(8 * vmax) - 0.5
    best = (np.inf, None)
    for r in np.linspace(0, R, 801):
        v = (r + 0.5) ** 2 / 8
        for th in np.linspace(0, 1, 41):
            W = th * r * r
            col = X.col_matrix(eqs, [(W, v)], S)[:, 0] @ y + (X.col_matrix(ineqs, [(W, v)], S)[:, 0] @ z if ni else 0)
            if col < best[0]:
                best = (col, (r, th, v, W))
    print(f"S = {S}: min Psi on fine grid {best[0]:.4e} at r, theta, v, W = {best[1]}")
