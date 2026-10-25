#!/usr/bin/env python3
"""src/cs/e3b/advchk/genus.py -- adversarial check of the E3b arithmetic-genus claim.
Korselt fibre for fixed (A, r_2..r_k): K_j (j = 2..k) is multi-affine in (c, l_i : i != j), degree 0 in l_j.
Closure in (P^1)^k (coords c, l_2..l_k): K_j has class D_j = h_1 + sum_{i>=2, i != j} h_i.
Adjunction for the complete intersection curve C = D_2 ... D_k:  K_C = (sum_j D_j - 2 sum_i h_i)|_C,
2 p_a - 2 = (sum_j D_j - 2 sum h) . D_2 ... D_k  in the ring Z[h_1..h_k]/(h_i^2)."""
import itertools, sys
from collections import defaultdict

def mul(P, Q):
    out = defaultdict(int)
    for a, ca in P.items():
        for b, cb in Q.items():
            if set(a) & set(b):
                continue          # h_i^2 = 0
            out[tuple(sorted(a + b))] += ca * cb
    return dict(out)

def genus(k):
    D = []
    for j in range(2, k + 1):
        D.append({(1,): 1, **{(i,): 1 for i in range(2, k + 1) if i != j}})
    Kc = defaultdict(int)
    for d in D:
        for m, cf in d.items():
            Kc[m] += cf
    for i in range(1, k + 1):
        Kc[(i,)] -= 2
    P = {(): 1}
    for d in D:
        P = mul(P, d)
    top = tuple(range(1, k + 1))
    deg_curve_classes = {i: mul(P, {(i,): 1}).get(top, 0) for i in range(1, k + 1)}
    P = mul(P, dict(Kc))
    degK = P.get(top, 0)
    return degK // 2 + 1, degK, deg_curve_classes, dict(Kc)

for k in (3, 5, 7):
    g, degK, dc, Kc = genus(k)
    print(f"k={k}: K_C class {Kc}  deg K_C = {degK}  arithmetic genus = {g}  "
          f"(degree of projection to each P^1 factor: {dc})")
