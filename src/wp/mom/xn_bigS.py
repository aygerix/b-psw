#!/usr/bin/env python3
"""xn line: the TRUE k = 9 problem at large S = sum v (some members near the linear member 2t + 1).

Member 1 is forced to have v_1 >= V0 (v_1 = V0 + exp(rho_1) ...); all members in the exact real domain.
Step 1: solve the four boundary equations (family + reflection) by least squares.
Step 2: maximise min(d2 fam, d2 refl) (d2 = 2^(k-1) T''(-1)) subject to them (SLSQP with slack).
A family needs d2 fam, d2 refl > 0.  COMPUTED only.
Usage: xn_bigS.py k lin V0 nstarts seed
"""
import sys
import numpy as np
from scipy.optimize import minimize, least_squares
from xn_core import vpoly_s

k, lin, V0, ns, seed = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
m = k - lin
rng = np.random.default_rng(seed)
sc = 2.0 ** (k - 1)


def unpack(p):
    r = np.exp(p[:m])
    v = (r + 0.5) ** 2 / 8
    v[0] += V0
    rr = np.sqrt(8 * v) - 0.5
    u = 0.5 + rr * np.tanh(p[m:])
    return u, v


def Tder(u, v):
    T = np.array(vpoly_s([float(x) for x in u], [float(x) for x in v], lin, k))
    n = np.arange(k)
    sg = (-1.0) ** n
    return np.dot(T, sg), np.dot(n[1:] * T[1:], sg[:-1]), np.dot((n * (n - 1))[2:] * T[2:], sg[2:])


def cons(p):
    u, v = unpack(p)
    a, b = Tder(u, v), Tder(1 - u, v)
    return np.array([sc * (a[0] - 2.0 ** (1 - k)), sc * a[1], sc * (b[0] - 2.0 ** (1 - k)), sc * b[1]])


def d2s(p):
    u, v = unpack(p)
    return np.array([sc * Tder(u, v)[2], sc * Tder(1 - u, v)[2]])


best = -np.inf
for it in range(ns):
    p0 = np.concatenate([rng.normal(-1, 1.5, m), rng.normal(0, 1.0, m)])
    with np.errstate(all="ignore"):
        try:
            s0 = least_squares(cons, p0, xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=3000)
        except Exception:
            continue
    if not np.all(np.isfinite(s0.fun)) or np.abs(s0.fun).max() > 1e-8:
        print(f"start {it}: boundary system not solved (maxres {np.abs(s0.fun).max():.2e})", flush=True)
        continue
    x0 = np.append(s0.x, d2s(s0.x).min())
    cs = [{"type": "eq", "fun": lambda q: cons(q[:-1])}, {"type": "ineq", "fun": lambda q: d2s(q[:-1]) - q[-1]}]
    with np.errstate(all="ignore"):
        try:
            r = minimize(lambda q: -q[-1], x0, constraints=cs, method="SLSQP", options={"maxiter": 1500, "ftol": 1e-12})
        except Exception:
            continue
    p = r.x[:-1]
    cv = np.abs(cons(p)).max()
    u, v = unpack(p)
    d = d2s(p)
    mg = (np.sqrt(8 * v) - (np.abs(u - 0.5) + 0.5)).min()
    print(f"start {it}: min d2 {d.min():.6g} (fam {d[0]:.6g}, refl {d[1]:.6g}) viol {cv:.1e} margin {mg:.2e} S {v.sum():.4f} "
          f"v sorted {np.array2string(np.sort(v)[::-1][:3], precision=4)}", flush=True)
    if cv < 1e-7 and d.min() > best:
        best = d.min()
print(f"k={k} lin={lin} V0={V0}: best min d2 = {best:.6g}")
