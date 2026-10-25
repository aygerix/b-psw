#!/usr/bin/env python3
"""xn line: maximise the second-order quantity d2 = 2^(k-1) T''(-1) subject to the boundary equations.

A real cubic Williams family has T(-1) = 2^(1-k), T'(-1) = 0 AND T''(-1) = 2^(-k) e_2(beta) > 0
(T(s) = s^(k-1) Vpoly(1/s); beta_i = h_i(-1) > 0), and the same for the reflection (beta -> gamma = g_i(0)).
mode 'fam'  : constraints E_0 = E_1 = 0 (family only), objective d2(family)
mode 'sum'  : constraints E0 + E0r = 0, E1 + E1r = 0 only (reflection-even part), objective d2 fam + d2 refl
mode 'both' : constraints E_0 = E_1 = 0 for family and reflection, objective min(d2 fam, d2 refl) (via slack)
If the maximum is < 0 then boundary equations + reality => contradiction with T''(-1) > 0.  COMPUTED only.
Usage: xn_maxd2.py k lin mode nstarts seed
"""
import sys
import numpy as np
from scipy.optimize import minimize, least_squares
from xn_core import vpoly_s
from xn_feas import unpack, resid, margins

k, lin, mode, ns, seed = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
m = k - lin
rng = np.random.default_rng(seed)
sc = 2.0 ** (k - 1)


def Tder(u, v):
    T = np.array(vpoly_s([float(x) for x in u], [float(x) for x in v], lin, k))
    n = np.arange(k)
    sg = (-1.0) ** n
    T0 = np.dot(T, sg)
    T1 = np.dot(n[1:] * T[1:], sg[:-1])
    T2 = np.dot((n * (n - 1))[2:] * T[2:], sg[2:])
    return T0, T1, T2


def parts(p):
    u, v = unpack(p)
    a = Tder(u, v)
    b = Tder(1 - u, v)
    return a, b


def cons(p):
    a, b = parts(p)
    c = [sc * (a[0] - 2.0 ** (1 - k)), sc * a[1]]
    if mode == "both":
        c += [sc * (b[0] - 2.0 ** (1 - k)), sc * b[1]]
    if mode == "sum":       # only the reflection-even combinations E0 + E0r, E1 + E1r
        c = [sc * (a[0] + b[0] - 2.0 ** (2 - k)), sc * (a[1] + b[1])]
    return np.array(c)


best = -np.inf
for it in range(ns):
    # start from a solution of the boundary system (least squares), then maximise
    p0 = np.concatenate([rng.normal(0, 1.5, m), rng.normal(0, 1.0, m)])
    with np.errstate(all="ignore"):
        try:
            s0 = least_squares(lambda p: cons(p), p0, xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=3000)
        except Exception:
            continue
    if not np.all(np.isfinite(s0.fun)) or np.abs(s0.fun).max() > 1e-8:
        continue
    if mode in ("fam", "sum"):
        x0 = s0.x
        obj = (lambda p: -sc * parts(p)[0][2]) if mode == "fam" else (lambda p: -sc * (parts(p)[0][2] + parts(p)[1][2]))
        cs = [{"type": "eq", "fun": cons}]
    else:
        a, b = parts(s0.x)
        x0 = np.append(s0.x, min(sc * a[2], sc * b[2]))
        obj = lambda q: -q[-1]
        cs = [{"type": "eq", "fun": lambda q: cons(q[:-1])},
              {"type": "ineq", "fun": lambda q: np.array([sc * parts(q[:-1])[0][2] - q[-1], sc * parts(q[:-1])[1][2] - q[-1]])}]
    with np.errstate(all="ignore"):
        try:
            r = minimize(obj, x0, constraints=cs, method="SLSQP", options={"maxiter": 2000, "ftol": 1e-12})
        except Exception:
            continue
    p = r.x if mode in ("fam", "sum") else r.x[:-1]
    cv = np.abs(cons(p)).max()
    a, b = parts(p)
    u, v = unpack(p)
    mg, dm = margins(u, v)
    val = sc * a[2] if mode == "fam" else (sc * (a[2] + b[2]) if mode == "sum" else min(sc * a[2], sc * b[2]))
    print(f"start {it}: d2 fam {sc*a[2]:.6g} refl {sc*b[2]:.6g} | constraint viol {cv:.1e} margin {mg:.2e} dmin {dm:.2e} S {v.sum():.5f} maxv {v.max():.4g}", flush=True)
    if cv < 1e-7:
        if val > best:
            best = val
            bu, bv = u.copy(), v.copy()
print(f"k={k} lin={lin} mode={mode}: best feasible objective {best:.6g}")
if np.isfinite(best):
    o = np.argsort(bv)
    print("   u =", np.array2string(bu[o], precision=8, max_line_width=300))
    print("   v =", np.array2string(bv[o], precision=8, max_line_width=300))
