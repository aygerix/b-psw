#!/usr/bin/env python3
"""src/cs/e3b/inspect_pt.py -- E3b: print, for a solution (c, l_2..l_k, r_2..r_k) of a stratum file, every root of
N_j with the relative Lucas residual |Phi_j|/|g_j|^(k-1), the Korselt residuals, and Lucas-point coincidences.
Usage: inspect_pt.py file k [filter: 'd>=1' style python expression on (d, gd)]"""
import sys, numpy as np
fn, k = sys.argv[1], int(sys.argv[2]); flt = sys.argv[3] if len(sys.argv) > 3 else "True"
seen = []
for line in open(fn):
    d, gd, clr, pts = eval(line)
    if not eval(flt):
        continue
    clr = np.array(clr)
    key = np.round(sorted(np.abs(clr)), 4).tolist()
    if key in seen:
        continue
    seen.append(key)
    c = clr[0]; lam = clr[1:k]; rr = clr[k:]
    g = lambda i, t: c if i == 0 else lam[i - 1] * (t - rr[i - 1])
    print("==== d", d, "gd", gd, "c,l,r =", np.round(clr, 6).tolist())
    for j in range(1, k):
        a = rr[j - 1]**2 - 1
        y = [g(i, rr[j - 1]) for i in range(k) if i != j]
        e = np.poly(-np.array(y))   # [1, e_1, ..., e_{k-1}]
        K = sum(a**(mm - 1) * e[mm] for mm in range(1, k))
        print(f"   Korselt {j+1}: |K| = {abs(K):.2e}   A(r) = {a:.5f}")
    for j in range(k):
        Ncoef = [c, 0, 2 - c] if j == 0 else [lam[j-1], -lam[j-1]*rr[j-1], -lam[j-1], lam[j-1]*rr[j-1] + 2]
        rts = np.roots(Ncoef)
        res = []
        for t in rts:
            val = np.prod([g(j, t) - 2 * g(i, t) for i in range(k) if i != j]) - g(j, t)**(k - 1)
            res.append((np.round(t, 5), f"{abs(val)/abs(g(j,t))**(k-1):.1e}"))
        print(f"   member {j+1}: roots of N_j and Lucas residuals: {res}")
