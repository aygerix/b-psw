#!/usr/bin/env python3
"""zsolve on adversarial F_5^r pools of size exactly (r+1)^5: must succeed (deterministic-bound claim).
Kinds: lowrank (span of r-1 vectors), repeats (3 distinct vectors), zeros (5% zero vectors), unit (e_j multiples),
twoval (each vector has entries in {0, c} for one c).  Usage: python3 adv_zsolve.py rmax trials seed"""
import sys, time
sys.path.insert(0, '..')
import numpy as np
import zsolve as z

def main():
    rmax, trials, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    g = np.random.default_rng(seed)
    for r in range(2, rmax + 1, 2):
        B = (r + 1) ** 5
        for kind in ("lowrank", "repeats", "zeros", "unit", "twoval"):
            used, t0 = [], time.time()
            for _ in range(trials):
                if kind == "lowrank":
                    Bas = g.integers(0, 5, size=(max(1, r - 1), r)); X = (g.integers(0, 5, size=(B, Bas.shape[0])) @ Bas) % 5
                elif kind == "repeats":
                    V = g.integers(1, 5, size=(3, r)); X = V[g.integers(0, 3, size=B)]
                elif kind == "zeros":
                    X = g.integers(0, 5, size=(B, r)); X[g.random(B) < 0.05] = 0
                elif kind == "unit":
                    X = np.zeros((B, r), dtype=np.int64); X[np.arange(B), g.integers(0, r, size=B)] = g.integers(1, 5, size=B)
                else:
                    X = (g.integers(0, 2, size=(B, r)) * g.integers(1, 5, size=(B, 1))) % 5
                XL = [tuple(int(a) for a in row) for row in X]
                S = z.Solver(z.Group([5] * r), XL); I = S.zero_sum()
                assert I is not None, f"FAILED within bound: r={r} {kind}"
                assert len(set(I)) == len(I) and not (X[np.array(I)].sum(axis=0) % 5).any()
                used.append(S.pulled)
            print(f"r={r:2d} {kind:8s} bound={B:8d} pulled max={max(used):8d} mean={sum(used)/len(used):10.1f} [{time.time()-t0:.1f}s]", flush=True)

if __name__ == "__main__":
    main()
