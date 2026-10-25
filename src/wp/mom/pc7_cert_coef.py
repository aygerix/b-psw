"""Coefficientwise certificates: D_alpha = kap + sum_m A_m p_m + sum_n c_n h_n + lam3 e3, where p_m, h_n
(h_n = sum_{i<j} sum_{a+b=n} w_i^a w_j^b) and e3 are > 0 for w > 0.  If some alpha makes kap_alpha > 0 and every A, c, lam3 >= 0,
the D-system has no solution with w > 0 at that shape.  Map of where this works, and whether D1 ALONE works.
Usage: python3 pc7_cert_coef.py step"""
import sys, numpy as np
from scipy.optimize import linprog
from pc7_cert_dlib import load, LAM3
at = load(); step = float(sys.argv[1])
vals = np.round(np.arange(-0.99, 0.991, step), 4)
ok = bad = d1 = 0; badl = []; d1l = []
for a in vals:
    for b in vals:
        if b < a or a + b < -1e-9: continue
        K, A, C = at(a, b)
        M = np.vstack([A.T, C.T, LAM3[None, :]])            # rows: coefficient functionals (15) + lam3, each must be >= 0
        Mn = M / np.maximum(np.abs(M).max(axis=1, keepdims=True), 1e-300)
        # maximise t: Mn alpha >= 0, kap alpha >= t, |alpha|_inf <= 1
        kn = K / np.abs(K).max()
        Aub = np.vstack([np.hstack([-Mn, np.zeros((len(Mn), 1))]), np.concatenate([-kn, [1]])[None, :]])
        r = linprog(np.concatenate([np.zeros(8), [-1]]), A_ub=Aub, b_ub=np.zeros(len(Aub)), bounds=[(-1, 1)] * 8 + [(None, None)], method='highs')
        good = r.status == 0 and r.x[-1] > 1e-9
        ok += good; bad += (not good)
        if not good: badl.append((a, b))
        if K[0] > 0 and np.all(M[:, 0] >= 0): d1 += 1; d1l.append((a, b))
print(f"step {step}: shapes {ok + bad}; coefficientwise certificate exists at {ok}; fails at {bad}; D1 alone coefficientwise-positive at {d1}")
print("fails (first 60):", badl[:60])
print("D1-alone region: r1 range", min(x for x, y in d1l) if d1l else None, "...", " sample:", d1l[:: max(1, len(d1l) // 30)])
