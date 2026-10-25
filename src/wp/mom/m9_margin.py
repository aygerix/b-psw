#!/usr/bin/env python3
"""m9 line: feasible points of the four boundary equations with LARGE interior margin (k = 9, 11).

Mode 'gen': m cubic members (u_i, v_i) free; 4 equations.
Mode 'sym': reflection-invariant configuration: p pairs (1/2 +- w_j, v_j), w_j > 0, plus (if m odd) a centre
member (1/2, v_c).  Its reflection is the same multiset, so the 4 equations reduce to the 2 family equations
(all 4 residuals are still evaluated and reported).
Phase 1: least squares in the exact domain parametrisation r = exp(rho), v = (r + 1/2)^2/8, w = r tanh(z)
(every point is admissible).  Phase 2 (from a phase-1 solution, feasible start): SLSQP maximising tau subject to
equations = 0, sqrt(8 v_i) - 1/2 - |w_i| >= tau, pairwise member distance >= tau/4 (and w_j >= tau/4 in 'sym').
A point is reported OK only if tau > 0, residual < 1e-9 and all inequalities hold.
Usage: m9_margin.py k lin mode nstarts seed
"""
import sys
import numpy as np
from scipy.optimize import minimize, least_squares
sys.path.insert(0, '.')
from m9_feas import Tvals


def nfree(m, mode):
    return m if mode == 'gen' else m // 2  # number of (v, w) blocks; plus centre in sym mode if m odd


def members(x, m, mode):
    """x = (v-block, w-block[, v_c]) ; returns u, v of all m cubic members"""
    if mode == 'gen':
        v, w = x[:m], x[m:2 * m]
        return 0.5 + w, v
    p = m // 2
    v, w = x[:p], x[p:2 * p]
    u = np.concatenate([0.5 + w, 0.5 - w])
    vv = np.concatenate([v, v])
    if m % 2:
        u = np.append(u, 0.5)
        vv = np.append(vv, x[2 * p])
    return u, vv


def from_param(q, m, mode):
    p = nfree(m, mode)
    r = np.exp(q[:p])
    v = (r + 0.5) ** 2 / 8
    w = r * np.tanh(q[p:2 * p]) if mode == 'gen' else r * (0.5 + 0.5 * np.tanh(q[p:2 * p]))
    out = [v, w]
    if mode == 'sym' and m % 2:
        out.append(np.array([(np.exp(q[2 * p]) + 0.5) ** 2 / 8]))
    return np.concatenate(out)


def eqres(x, k, lin, m, mode):
    u, v = members(x, m, mode)
    c0 = 2.0 ** (1 - k)
    a0, a1 = Tvals(u, v, k, lin)
    b0, b1 = Tvals(1 - u, v, k, lin)
    return np.array([(a0 - c0) / c0, a1 / c0, (b0 - c0) / c0, b1 / c0])


def ineqs(x, m, mode):
    u, v = members(x, m, mode)
    mg = np.sqrt(8 * np.maximum(v, 0.0)) - 0.5 - np.abs(u - 0.5)
    d = np.array([np.hypot(u[i] - u[j], v[i] - v[j]) for i in range(m) for j in range(i + 1, m)])
    ws = x[m // 2:2 * (m // 2)] if mode == 'sym' else np.zeros(0)
    return mg, d, ws


def main():
    k, lin, mode, ns, seed = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
    m = k - lin
    p = nfree(m, mode)
    nq = 2 * p + (1 if (mode == 'sym' and m % 2) else 0)
    rng = np.random.default_rng(seed)
    best = None
    neq = 2 if mode == 'sym' else 4
    for it in range(ns):
        q0 = np.concatenate([rng.normal(-0.5, 1.2, p), rng.normal(0, 1.0, nq - p)])
        try:
            s1 = least_squares(lambda q: eqres(from_param(q, m, mode), k, lin, m, mode)[:neq], q0,
                               xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
        except Exception as e:
            print(f"start {it}: phase1 exception {e!r}", flush=True)
            continue
        x1 = from_param(s1.x, m, mode)
        if np.abs(eqres(x1, k, lin, m, mode)).max() > 1e-8:
            print(f"start {it}: phase1 fail {np.abs(s1.fun).max():.2e}", flush=True)
            continue
        mg, d, ws = ineqs(x1, m, mode)
        tau0 = min(mg.min(), 4 * d.min(), 4 * ws.min() if len(ws) else 1e9)

        def cin(z):
            mg, d, ws = ineqs(z[:-1], m, mode)
            return np.concatenate([mg - z[-1], d - z[-1] / 4, ws - z[-1] / 4])
        cons = [{'type': 'eq', 'fun': lambda z: eqres(z[:-1], k, lin, m, mode)[:neq]},
                {'type': 'ineq', 'fun': cin}]
        try:
            res = minimize(lambda z: -z[-1], np.append(x1, tau0), method='SLSQP', constraints=cons,
                           options={'maxiter': 3000, 'ftol': 1e-15})
            z = res.x
        except Exception as e:
            print(f"start {it}: phase2 exception {e!r}", flush=True)
            z = np.append(x1, tau0)
        er = np.abs(eqres(z[:-1], k, lin, m, mode)).max()
        mg, d, ws = ineqs(z[:-1], m, mode)
        tau = min(mg.min(), 4 * d.min(), 4 * ws.min() if len(ws) else 1e9)
        ok = er < 1e-9 and tau > 0
        if not ok:   # fall back to the phase-1 point
            z, tau = np.append(x1, tau0), tau0
            ok = tau0 > 0
        print(f"start {it}: phase1 tau {tau0:.5f} -> tau {tau:.5f} eqres {er:.2e} {'OK' if ok else ''}", flush=True)
        if ok and (best is None or tau > best[0]):
            best = (tau, z.copy())
    if best is not None:
        tau, z = best
        u, v = members(z[:-1], m, mode)
        np.set_printoptions(precision=16, linewidth=250)
        print("BEST tau =", tau)
        print("x =", repr(z[:-1]))
        print("u =", repr(u))
        print("v =", repr(v))
        print("residuals =", eqres(z[:-1], k, lin, m, mode))
        mg, d, ws = ineqs(z[:-1], m, mode)
        print("margins sqrt(8v)-1/2-|w| =", mg, " min dist", d.min())


if __name__ == "__main__":
    main()
