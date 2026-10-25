#!/usr/bin/env python3
"""Scaling of the exponent-5 solver (norm trick + two doublings + LA): elements pulled for one 0/1 zero-sum in
(Z/5)^r versus the deterministic bound (r+1)^5; also the ternary (QR) stage alone versus (r+1)^2."""
import sys, random, time, math
import numpy as np
sys.path.insert(0, '..')
import zsolve as z

def ternary_only(r, rng):
    X = [np.array([rng.randrange(5) for _ in range(r)]) for _ in range((r + 1) ** 2 + 10)]
    it = iter(range(len(X))); cnt = [0]
    def gen():
        i = next(it, None)
        if i is None: return None
        cnt[0] = i + 1
        return z.Elt(X[i], {1: ((i,), (0,))})
    rel = z.QRRel(z.Stream(gen), 5, r).pull()
    tot = sum(c * X[e.reals[1][0][0]] for e, c in rel.items()) % 5
    assert not tot.any() and all(c in (1, 4) for c in rel.values())
    return cnt[0], len(rel)

rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 11)
for r in range(2, int(sys.argv[1]) + 1, 2):
    tr = 5 if r <= 10 else 2
    t0 = time.time(); pulls = []; tern = []; sizes = []
    for _ in range(tr):
        G = z.Group([5] * r); B = (r + 1) ** 5
        X = [tuple(rng.randrange(5) for _ in range(r)) for _ in range(min(B, 2_000_000))]
        S = z.Solver(G, X); I = S.zero_sum(); assert I
        assert len(set(I)) == len(I) and all(0 <= i < len(X) for i in I), 'not a 0/1 subset'
        assert not (np.array([X[i] for i in I]).sum(axis=0) % 5).any(), 'numpy recheck failed'
        tot = G.zero()
        for i in I: tot = G.add(tot, X[i])
        assert tot == G.zero()
        sizes.append(len(I)); pulls.append(S.pulled); tern.append(ternary_only(r, rng)[0])
    m = sum(pulls) / tr
    print(f"r={r:3d} ratio={m/(r+1)**5:.3f} min={min(pulls)} max={max(pulls)} bound={(r+1)**5:10d} pulled mean={m:10.1f} log_(r+1)={math.log(m)/math.log(r+1):5.2f} "
          f"|I|={sizes} pulls={pulls} | ternary stage: pulled mean {sum(tern)/tr:8.1f} vs (r+1)^2={(r+1)**2}  [{time.time()-t0:.0f}s]", flush=True)
