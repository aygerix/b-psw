#!/usr/bin/env python3
import os
"""Run the ORIGINAL smallexp/zsolve.py solver on adversarial (Z/3)^r inputs of length exactly (r+1)^2;
verify outputs.  Also print sq_table(3) to confirm t = a - 2b with a = support, b = {t = 2}."""
import sys, random
sys.path.insert(0, os.path.expanduser("~/b-psw/src/cs/smallexp"))
sys.path.insert(0, os.path.expanduser("~/b-psw/src/cs/audit_l3"))
import zsolve as z
import l3_indep as ind
print("sq_table(3) =", z.sq_table(3))
rng = random.Random(99)
for r in (1, 2, 3, 5, 8, 13, 20):
    G = z.Group([3] * r); n = (r + 1) ** 2
    for kind in ind.KINDS:
        for _ in range(100 if r <= 8 else 20):
            X = ind.gen(kind, r, n, rng)
            S = z.Solver(G, X)
            I = S.zero_sum()
            assert I is not None, (r, kind)
            assert S.pulled <= n
            ind.check(X, r, I)
    print(f"r={r}: original solver OK on all families within (r+1)^2 = {n}", flush=True)
