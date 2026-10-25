#!/usr/bin/env python3
"""m9 line: k = 11 (eleven cubics) in the reflection-symmetric slice, seeded from the k = 9 symmetric solution.

Slice: p pairs (1/2 +- w_j, v_j) + one self-dual member (1/2, v_c); reflection-invariant, so only the two family
boundary equations remain.  Start: the k = 9 symmetric solution (4 pairs + centre) plus one new random pair;
least squares on the k = 11 equations (exact domain parametrisation).  Optionally also a general refinement.
Usage: m9_symseed.py ntrials seed
"""
import sys
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, '.')
from m9_feas import Tvals

W9 = np.array([1.6286121386777186, 0.0286312273614805, 0.1347922858160165, 0.1143013292143569])
V9 = np.array([0.605156204392249, 0.053680237637796, 0.0620060685662747, 0.059648427249592])
VC9 = 1.6107368111675027


def build(q, p):
    r = np.exp(q[:p])
    v = (r + .5) ** 2 / 8
    w = r * (1 + np.tanh(q[p:2 * p])) / 2
    vc = (np.exp(q[2 * p]) + .5) ** 2 / 8
    u = np.concatenate([.5 + w, .5 - w, [.5]])
    return u, np.concatenate([v, v, [vc]])


def tobuild(w, v, vc):
    r = np.sqrt(8 * v) - .5
    z = np.arctanh(np.clip(2 * w / r - 1, -0.999999, 0.999999))
    return np.concatenate([np.log(r), z, [np.log(np.sqrt(8 * vc) - .5)]])


def resid(q, p, k):
    u, v = build(q, p)
    c0 = 2.0 ** (1 - k)
    a0, a1 = Tvals(u, v, k, 0)
    b0, b1 = Tvals(1 - u, v, k, 0)
    return np.array([(a0 - c0) / c0, a1 / c0, (b0 - c0) / c0, b1 / c0])


def main():
    nt, seed = map(int, sys.argv[1:3])
    rng = np.random.default_rng(seed)
    p = 5
    nok = 0
    for it in range(nt):
        v0 = 1 / 32 + np.exp(rng.normal(-3, 1.5))
        w0 = (np.sqrt(8 * v0) - .5) * rng.uniform(0.05, 0.95)
        w = np.append(W9, w0) * np.exp(rng.normal(0, 0.05, 5))
        v = np.append(V9, v0) * np.exp(rng.normal(0, 0.05, 5))
        v = np.maximum(v, (w + .5) ** 2 / 8 * 1.0001)
        q0 = tobuild(w, v, VC9 * np.exp(rng.normal(0, 0.1)))
        try:
            s = least_squares(lambda q: resid(q, p, 11)[:2], q0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
        except Exception:
            continue
        u, vv = build(s.x, p)
        f = np.abs(resid(s.x, p, 11)).max()
        mg = (np.sqrt(8 * vv) - .5 - np.abs(u - .5)).min()
        dm = min(np.hypot(u[i] - u[j], vv[i] - vv[j]) for i in range(len(u)) for j in range(i + 1, len(u)))
        ok = f < 1e-8 and mg > 1e-6 and dm > 1e-4
        print(f"trial {it}: maxres {f:.2e} margin {mg:.2e} dmin {dm:.2e} {'SOL' if ok else 'fail'}", flush=True)
        if ok:
            nok += 1
            np.set_printoptions(precision=16, linewidth=250)
            print("   u =", repr(u))
            print("   v =", repr(vv))
    print("solutions:", nok, "of", nt)


if __name__ == "__main__":
    main()
