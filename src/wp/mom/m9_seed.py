#!/usr/bin/env python3
"""m9 line: k = 11 boundary-equation feasibility, seeded from k = 9 solutions.

Loop: (1) find a k = 9 solution (same lin) by least squares from a random start (exact domain parametrisation);
(2) append two members (random type: corner-ish small v, pair-like, self-dual, or random) and solve the k = 11
equations by least squares starting from that configuration.  Prints k = 11 solutions with margin and distance.
Usage: m9_seed.py lin ntrials seed
"""
import sys
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, '.')
from m9_feas import Tvals


def unpack(q, m):
    r = np.exp(q[:m])
    return 0.5 + r * np.tanh(q[m:]), (r + 0.5) ** 2 / 8


def pack(u, v):
    r = np.sqrt(8 * v) - 0.5
    z = np.arctanh(np.clip((u - 0.5) / r, -0.999999, 0.999999))
    return np.concatenate([np.log(r), z])


def resid(q, k, lin):
    m = k - lin
    u, v = unpack(q, m)
    c0 = 2.0 ** (1 - k)
    a0, a1 = Tvals(u, v, k, lin)
    b0, b1 = Tvals(1 - u, v, k, lin)
    return np.array([(a0 - c0) / c0, a1 / c0, (b0 - c0) / c0, b1 / c0])


def solve(q0, k, lin):
    s = least_squares(resid, q0, args=(k, lin), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
    return s, np.abs(s.fun).max()


def main():
    lin, nt, seed = map(int, sys.argv[1:4])
    rng = np.random.default_rng(seed)
    n9 = n11 = 0
    for it in range(nt):
        m9 = 9 - lin
        q0 = np.concatenate([rng.normal(-0.5, 1.2, m9), rng.normal(0, 1.0, m9)])
        try:
            s, f = solve(q0, 9, lin)
        except Exception:
            continue
        if f > 1e-8:
            continue
        n9 += 1
        u9, v9 = unpack(s.x, m9)
        for rep in range(4):
            typ = rng.integers(4)
            if typ == 0:   # two corner-ish members
                vv = 1 / 32 + rng.uniform(0.005, 0.05, 2)
                uu = 0.5 + (np.sqrt(8 * vv) - 0.5) * rng.uniform(-0.8, 0.8, 2)
            elif typ == 1:  # a reflected pair
                v0 = (np.exp(rng.normal(-0.5, 1)) + .5) ** 2 / 8
                w0 = (np.sqrt(8 * v0) - .5) * rng.uniform(0.1, 0.9)
                vv, uu = np.array([v0, v0]), np.array([.5 + w0, .5 - w0])
            elif typ == 2:  # two self-dual members
                vv = (np.exp(rng.normal(-0.5, 1.2, 2)) + .5) ** 2 / 8
                uu = np.array([.5, .5])
            else:
                r = np.exp(rng.normal(-0.5, 1.2, 2))
                vv, uu = (r + .5) ** 2 / 8, .5 + r * rng.uniform(-.9, .9, 2)
            u, v = np.concatenate([u9, uu]), np.concatenate([v9, vv])
            try:
                s2, f2 = solve(pack(u, v), 11, lin)
            except Exception:
                continue
            u2, v2 = unpack(s2.x, 11 - lin)
            mg = (np.sqrt(8 * v2) - .5 - np.abs(u2 - .5)).min()
            dm = min(np.hypot(u2[i] - u2[j], v2[i] - v2[j]) for i in range(len(u2)) for j in range(i + 1, len(u2)))
            if f2 < 1e-8:
                n11 += 1
                np.set_printoptions(precision=16, linewidth=250)
                print(f"trial {it}.{rep} type {typ}: k=11 SOL maxres {f2:.1e} margin {mg:.5f} dmin {dm:.5f}", flush=True)
                print("   u =", repr(u2))
                print("   v =", repr(v2))
            else:
                print(f"trial {it}.{rep} type {typ}: k=11 fail {f2:.2e} (margin {mg:.1e}, dmin {dm:.1e})", flush=True)
    print("k=9 seeds:", n9, " k=11 solutions:", n11)


if __name__ == "__main__":
    main()
