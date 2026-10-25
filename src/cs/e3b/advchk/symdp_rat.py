#!/usr/bin/env python3
import os
"""src/cs/e3b/advchk/symdp_rat.py -- field test for the REAL Sum m = 6 families on the symmetric double-pencil component
(out/symdp6.txt from symdp_nodes.py).  Polish to dps digits (mpmath Newton on the square node system), then test whether the
normalised invariants c, l_2 r, l_3 r, r^2 are rational (mpmath.findpoly degree 1, |coeff| <= 1e15); a family over Q needs all
four rational (A = x^2 - 1 normalisation, Lemma-B-style argument).  Also re-counts Sum m at high precision.
Usage: symdp_rat.py file [dps]"""
import sys, numpy as np, mpmath as mp
sys.path.insert(0, os.path.expanduser('~/b-psw/src/cs/e3b/advchk'))
from symdp_nodes import SymDP
fn = sys.argv[1]; dps = int(sys.argv[2]) if len(sys.argv) > 2 else 60
mp.mp.dps = dps
models = {}
nreal = 0; nrat = 0; worst = mp.mpf(0)
for line in open(fn):
    st, gd, real, z = eval(line)
    if not real:
        continue
    nreal += 1
    S = models.setdefault(st, SymDP(st))
    def Fv(w):
        Z = np.array([w], dtype=object)
        return [mp.mpc(v) for v in S.full(Z)[0]]
    w = [mp.mpc(v) for v in z]
    n = len(w)
    for it in range(60):
        f = Fv(w); nf = mp.norm(mp.matrix(f))
        if nf < mp.mpf(10) ** (-dps + 8):
            break
        h = mp.mpf(10) ** (-dps // 2)
        J = mp.matrix(n, n)
        for i in range(n):
            wp = list(w); wp[i] += h; wm = list(w); wm[i] -= h
            fp = Fv(wp); fm = Fv(wm)
            for r_ in range(n):
                J[r_, i] = (fp[r_] - fm[r_]) / (2 * h)
        d = mp.lu_solve(J, -mp.matrix(f))
        w = [w[i] + d[i] for i in range(n)]
    worst = max(worst, nf)
    c, l2, l3, r = w[:4]
    inv = [mp.re(c), mp.re(l2 * r), mp.re(l3 * r), mp.re(r * r)]
    rat = [bool(mp.findpoly(v, 1, maxcoeff=10**15)) for v in inv]
    if all(rat):
        nrat += 1
        print("RATIONAL FAMILY:", st, gd, [mp.nstr(v, 30) for v in inv], flush=True)
    elif any(rat):
        print("partially rational invariants:", st, gd, rat, [mp.nstr(v, 20) for v in inv], flush=True)
print(f"real families tested: {nreal}; fully rational: {nrat}; worst polish residual {mp.nstr(worst, 3)}")
