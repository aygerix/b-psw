#!/usr/bin/env python3
"""xn line: at solutions of the boundary system (N = 2), evaluate the next necessary conditions.

For each N = 2 solution (family and reflection):
  d2V  := 2^(k-1) Vpoly''(-1) - (k-1)(k-2)   [must equal e_2(beta)/2 > 0, beta_i = 2 + a_i - c_i = h_i(-1)]
  e2b  := e_2(beta)/2 at the point
  minD := min over real t of D(t)/(t+1)^2, D = Vpoly - 2^(1-k) t^(k-1)   [must be >= 0 for a family]
  Vmin := min over real t of Vpoly(t)                                       [must be > 0]
Usage: xn_n2diag.py k lin nstarts seed
"""
import sys
import numpy as np
from scipy.optimize import least_squares
from xn_feas import unpack, resid, margins
from xn_core import vpoly_s

k, lin, ns, seed = map(int, sys.argv[1:5])
m = k - lin
rng = np.random.default_rng(seed)


def vpoly(u, v):
    T = vpoly_s([float(x) for x in u], [float(x) for x in v], lin, k)
    return np.array(T[::-1])   # numpy poly1d order reversed: coefficient of t^0 first


def info(u, v):
    c = vpoly(u, v)                     # c[j] = coeff of t^j
    P = np.polynomial.Polynomial(c)
    d2 = 2.0 ** (k - 1) * P.deriv(2)(-1.0) - (k - 1) * (k - 2)
    a = 1 / v
    cc = u / v
    beta = 2 + a - cc
    if lin:
        beta = np.append(beta, 2.0)
    e2 = (beta.sum() ** 2 - (beta ** 2).sum()) / 2
    D = P - np.polynomial.Polynomial([0] * (k - 1) + [2.0 ** (1 - k)])
    Q, R = divmod(D, np.polynomial.Polynomial([1, 2, 1]))
    ts = np.linspace(-6, 6, 24001)
    return d2, e2 / 2, Q(ts).min(), P(ts).min(), np.abs(R.coef).max()


cnt = 0
for it in range(ns):
    p0 = np.concatenate([rng.normal(0, 1.5, m), rng.normal(0, 1.0, m)])
    with np.errstate(all="ignore"):
        try:
            sol = least_squares(resid, p0, args=(k, lin, 2), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        except Exception:
            continue
    if not np.all(np.isfinite(sol.fun)) or np.abs(sol.fun).max() > 1e-9:
        continue
    u, v = unpack(sol.x)
    mg, dm = margins(u, v)
    if dm < 1e-6 or mg < 1e-8:
        continue
    cnt += 1
    f = info(u, v)
    r = info(1 - u, v)
    print(f"sol {cnt}: margin {mg:.2e} dmin {dm:.2e} | fam: d2V {f[0]:.6g} e2b/2 {f[1]:.6g} minD/(t+1)^2 {f[2]:.4g} minV {f[3]:.3g}"
          f" | refl: d2V {r[0]:.6g} e2g/2 {r[1]:.6g} minD {r[2]:.4g} minV {r[3]:.3g}", flush=True)
print("solutions:", cnt)
