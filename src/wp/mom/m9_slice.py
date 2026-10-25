#!/usr/bin/env python3
"""m9 line: numerical feasibility of the four boundary equations on structured slices (k = 9, 11).

Members: nA asymmetric (u, v) free; nP reflected pairs (1/2 +- w, v); nB self-dual members (1/2, v).
m = nA + 2 nP + nB cubic members (+ 2t+1 if lin).  Exact domain parametrisation (every point admissible):
r = exp(rho), v = (r + 1/2)^2 / 8, w = r tanh(z)   (pairs: w = r (1 + tanh z)/2 > 0; self-dual: w = 0).
Least squares on the 4 equations from random starts; prints solutions with min margin and min member distance.
Usage: m9_slice.py k lin nA nP nB nstarts seed
"""
import sys
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, '.')
from m9_feas import Tvals


def build(q, nA, nP, nB):
    i = 0
    ra = np.exp(q[i:i + nA]); i += nA
    za = q[i:i + nA]; i += nA
    rp = np.exp(q[i:i + nP]); i += nP
    zp = q[i:i + nP]; i += nP
    rb = np.exp(q[i:i + nB]); i += nB
    va, vp, vb = (ra + .5) ** 2 / 8, (rp + .5) ** 2 / 8, (rb + .5) ** 2 / 8
    wa, wp = ra * np.tanh(za), rp * (1 + np.tanh(zp)) / 2
    u = np.concatenate([.5 + wa, .5 + wp, .5 - wp, .5 * np.ones(nB)])
    v = np.concatenate([va, vp, vp, vb])
    return u, v


def resid(q, k, lin, nA, nP, nB):
    u, v = build(q, nA, nP, nB)
    c0 = 2.0 ** (1 - k)
    a0, a1 = Tvals(u, v, k, lin)
    b0, b1 = Tvals(1 - u, v, k, lin)
    return np.array([(a0 - c0) / c0, a1 / c0, (b0 - c0) / c0, b1 / c0])


def main():
    k, lin, nA, nP, nB, ns, seed = map(int, sys.argv[1:8])
    assert nA + 2 * nP + nB + lin == k
    rng = np.random.default_rng(seed)
    nq = 2 * nA + 2 * nP + nB
    nok = 0
    best = None
    for it in range(ns):
        q0 = rng.normal(0, 1.2, nq)
        q0[:nA] -= 0.5
        try:
            s = least_squares(resid, q0, args=(k, lin, nA, nP, nB), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
        except Exception:
            continue
        f = np.abs(s.fun).max()
        u, v = build(s.x, nA, nP, nB)
        mg = (np.sqrt(8 * v) - .5 - np.abs(u - .5)).min()
        dm = min(np.hypot(u[i] - u[j], v[i] - v[j]) for i in range(len(u)) for j in range(i + 1, len(u)))
        if f < 1e-9:
            nok += 1
            sc = min(mg, dm)
            print(f"start {it}: SOL maxres {f:.1e} margin {mg:.5f} dmin {dm:.5f}", flush=True)
            if best is None or sc > best[0]:
                best = (sc, u, v)
        else:
            print(f"start {it}: fail {f:.2e} (margin {mg:.1e}, dmin {dm:.1e})", flush=True)
    print("solutions:", nok, "of", ns)
    if best:
        np.set_printoptions(precision=14, linewidth=250)
        print("best min(margin, dmin) =", best[0])
        print("u =", repr(best[1]))
        print("v =", repr(best[2]))


if __name__ == "__main__":
    main()
