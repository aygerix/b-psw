#!/usr/bin/env python3
"""Audit of the layering claim for groups of exponent 2^a 3^b 5^c (uses zsolve.Solver; sums rechecked
independently with numpy, subsets checked to be 0/1, pulled <= deterministic bound, targets checked).
Also adversarial pools: elements confined to a proper subgroup, and pools rich in elements of order 5^c.
Usage: python3 layer_chk.py seed"""
import sys, time, random
sys.path.insert(0, '..')
import numpy as np
import zsolve as z

def check(mods, X, I, u=None):
    M = np.array(mods, dtype=np.int64)
    s = np.zeros(len(mods), dtype=np.int64) if u is None else np.array(u, dtype=np.int64)
    if I:
        s = s + np.array([X[i] for i in I], dtype=np.int64).sum(axis=0)
    assert len(set(I)) == len(I), "repeated index"
    assert not (s % M).any(), "sum not zero"

def run(mods, trials, g, kind):
    G = z.Group(mods); B = z.bound(G); N = int(min(B, 2_000_000))
    M = np.array(mods, dtype=np.int64)
    P, T, t0 = [], [], time.time()
    for _ in range(trials):
        if kind == "rand":
            X = [tuple(int(a) for a in row) for row in (g.integers(0, 10**9, size=(N, len(mods))) % M)]
        elif kind == "sub":      # elements in the subgroup 5*G (or 2*G if no 5): the shallow layers vanish
            X = [tuple(int(a) for a in row) for row in ((g.integers(0, 10**9, size=(N, len(mods))) * 5) % M)]
        else:                    # 'deep': each element has a single nonzero coordinate
            X = []
            for _ in range(N):
                v = [0] * len(mods); j = int(g.integers(len(mods))); v[j] = int(g.integers(1, mods[j])); X.append(tuple(v))
        S = z.Solver(G, X); I = S.zero_sum()
        assert I is not None, f"no zero-sum within bound {B} for {mods} {kind}"
        check(mods, X, I); assert S.pulled <= B; P.append(S.pulled)
        if kind == "rand":
            u = tuple(int(a) for a in g.integers(0, 10**9, size=len(mods)) % M)
            # a target can need up to (rho+1) more top-level elements per layer: give a 2B pool
            X2 = X + [tuple(int(a) for a in row) for row in (g.integers(0, 10**9, size=(N, len(mods))) % M)]
            S2 = z.Solver(G, X2); J = S2.target(u)
            assert J is not None, "target failed"
            check(mods, X2, J, u); T.append(S2.pulled)
    lay = [(l, k, len(a)) for l, k, a in G.layers()]
    print(f"{str(mods):22s} {kind:4s} layers={lay} bound={B:9d} pulled max={max(P):8d} mean={sum(P)/len(P):10.1f}"
          + (f" target mean={sum(T)/len(T):10.1f}" if T else "") + f" [{time.time()-t0:.1f}s]", flush=True)

def main():
    g = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
    cases = [([5] * 3, 5), ([25] * 2, 3), ([125], 10), ([125, 5], 3), ([50, 10], 3), ([75, 15], 3),
             ([60, 30], 3), ([150], 10), ([100, 20], 2), ([90, 45], 2), ([625], 5), ([5, 25], 2),
             ([2, 3, 5, 4, 9, 25], 1), ([20, 30, 45], 1)]
    for mods, tr in cases:
        for kind in ("rand", "sub", "deep"):
            run(mods, tr, g, kind)

if __name__ == "__main__":
    main()
