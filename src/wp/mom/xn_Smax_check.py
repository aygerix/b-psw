#!/usr/bin/env python3
"""xn line: numerical check of the identity behind the a priori bound (xn_Smax.gp) at solutions of the boundary
equations (which satisfy the even boundary equations used by the identity):
  lin 0: Re(-7/5) - 59/110 Phi = -2/55 S^2 - 788/125 S - 2/55 B0 + 39391023/1100000
  lin 1: Re(-7/5) - 27/50  Phi = -1/25 S^2 - 678/125 S - 1/25 B0 + 1294023/50000
Re(t) = (Ve(t) - 2^-8 t^8)/(t+1)^2, Ve = (V + V~)/2 computed directly from the Laurent series (xn_core.vpoly_s);
Phi = (D2 + D2r)/2.  Also reports which ingredient fails when S > Smax (Re(-7/5) < 0 or Phi >= PhiBar)."""
import sys
import numpy as np
from scipy.optimize import least_squares
from xn_feas import unpack, resid, margins
from xn_core import vpoly_s

lin, ns, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
k, m = 9, 9 - lin
rng = np.random.default_rng(seed)
lam, q2, q1, qB, q0 = ((59 / 110, -2 / 55, -788 / 125, -2 / 55, 39391023 / 1100000) if lin == 0 else
                       (27 / 50, -1 / 25, -678 / 125, -1 / 25, 1294023 / 50000))
PhiBar = 729 / 32 if lin == 0 else 585 / 32
cnt = 0
for it in range(ns):
    p0 = np.concatenate([rng.normal(0, 1.5, m), rng.normal(0, 1.0, m)])
    with np.errstate(all="ignore"):
        try:
            sol = least_squares(resid, p0, args=(k, lin, 2), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
        except Exception:
            continue
    if not np.all(np.isfinite(sol.fun)) or np.abs(sol.fun).max() > 1e-9:
        continue
    u, v = unpack(sol.x)
    mg, dm = margins(u, v)
    if mg <= 0 or dm < 1e-6:
        continue
    cnt += 1
    T = np.array(vpoly_s(list(u), list(v), lin, k))
    Tr = np.array(vpoly_s(list(1 - u), list(v), lin, k))
    Te = (T + Tr) / 2
    Ve = lambda t: sum(Te[n] * t ** (8 - n) for n in range(9))
    d2 = lambda TT: sum(TT[n] * (8 - n) * (7 - n) * (-1.0) ** n for n in range(7))
    Phi = (d2(T) + d2(Tr)) / 2 - 56 / 256
    t0 = -1.4
    Re = (Ve(t0) - t0 ** 8 / 256) / (t0 + 1) ** 2
    S, B0 = v.sum(), (v ** 2).sum()
    lhs = Re - lam * Phi
    rhs = q2 * S * S + q1 * S + qB * B0 + q0
    print(f"sol {cnt}: S {S:.5f}  identity residual {lhs - rhs:.3e} (scale {abs(rhs):.3e})   Re(-7/5) {Re:.4g}  Phi {Phi:.4g} (PhiBar {PhiBar:.4g})", flush=True)
