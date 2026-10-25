#!/usr/bin/env python3
"""Audit of zsolve.QRRel (the claim's code): (1) adversarial inputs, EXACTLY (r+1)^2 elements supplied
(tests the deterministic bound), several odd primes; (2) l = 5 ternary-stage consumption on uniform inputs,
more trials and larger r than scale5.py. Relations re-verified from the raw vectors."""
import sys, random, time, statistics
import numpy as np
sys.path.insert(0, '..')
import zsolve as z
from qr_indep import gens, issq

def run_code(X, l, r):
    it = iter(range(len(X))); mx = [0]
    def gen():
        i = next(it, None)
        if i is None: return None
        mx[0] = max(mx[0], i + 1)
        return z.Elt(np.array(X[i], dtype=np.int64), {1: ((i,), (0,))})
    rel = z.QRRel(z.Stream(gen), l, r).pull()
    out = {e.reals[1][0][0]: int(c) for e, c in rel.items()}
    assert len(out) == len(rel) and out
    assert all(c % l and issq(c, l) for c in out.values())
    tot = [sum(c * X[i][k] for i, c in out.items()) % l for k in range(r)]
    assert not any(tot)
    return mx[0], out

if __name__ == '__main__':
    mode = sys.argv[1]
    rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 3)
    if mode == 'adv':
        n = 0
        for l in [3, 5, 7, 11, 13, 17, 29, 31]:
            for r in [1, 2, 3, 5, 8, 13]:
                for _ in range(3):
                    for name, X in gens(l, r, rng):
                        try:
                            used, rel = run_code(X[:(r + 1) ** 2], l, r)
                        except z.Exhausted:
                            print("EXHAUSTED within (r+1)^2:", l, r, name); raise
                        assert used <= (r + 1) ** 2
                        n += 1
        print("adversarial code runs passed:", n)
    else:
        trials = int(sys.argv[3]) if len(sys.argv) > 3 else 20
        for r in [8, 12, 16, 20, 24, 32, 40]:
            t0 = time.time(); us = []; sz = []; pm = []
            for _ in range(trials):
                X = [[rng.randrange(5) for _ in range(r)] for _ in range((r + 1) ** 2)]
                used, rel = run_code(X, 5, r)
                us.append(used); sz.append(len(rel)); pm.append(sum(1 for c in rel.values() if c == 1))
            print(f"l=5 r={r:3d} (r+1)^2={(r+1)**2:5d} used mean {statistics.mean(us):8.1f} "
                  f"min {min(us)} max {max(us)} ratio {statistics.mean(us)/(r+1)**2:.3f} | |rel| mean "
                  f"{statistics.mean(sz):6.1f}, #(+1) mean {statistics.mean(pm):6.1f} [{time.time()-t0:.0f}s]", flush=True)
