"""Decomposed-positivity LP at fixed shapes: find alpha in R^8 with
   kap(alpha) = 1, a_alpha(w) >= 0 (w > 0), b_alpha(w, v) >= 0 (w, v > 0), lam3(alpha) >= 0,
which makes sum alpha_i E_i > 0 for every w in (0, inf)^7.  Usage: python3 pc7_cert_dp.py ngrid [mode]"""
import sys, time, numpy as np
from scipy.optimize import linprog
from pc7_cert_lib import load, at_shape, LAM3, NAMES
D = load()
WG = np.exp(np.linspace(np.log(1e-6), np.log(1e6), 241))
PW = np.exp(np.linspace(np.log(1e-4), np.log(1e4), 41))
WV = np.array([(a, b) for a in PW for b in PW if b >= a])
def rows(R1, R2):
    kap, Ac, Bc = at_shape(D, R1, R2)
    Ra = np.array([[np.polyval(Ac[i][::-1], x) for i in range(8)] for x in WG])
    Rb = np.array([[np.sum(Bc[i] * np.outer(a ** np.arange(10), b ** np.arange(10))) for i in range(8)] for a, b in WV])
    Ra /= np.abs(Ra).max(axis=1, keepdims=True); Rb /= np.maximum(np.abs(Rb).max(axis=1, keepdims=True), 1e-300)
    return kap, Ac, Ra, Rb
def dp(R1, R2, mode='kap1'):
    kap, Ac, Ra, Rb = rows(R1, R2)
    # leading / trailing coefficient constraints (asymptotics w -> inf, 0)
    Aub = -np.vstack([Ra, Rb, LAM3[None, :]])
    bub = np.zeros(len(Aub))
    if mode == 'kap1': Aeq, beq = kap[None, :], [1.0]
    res = linprog(np.zeros(8), A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=[(-1e4, 1e4)] * 8, method='highs')
    return res
ng = int(sys.argv[1]); mode = sys.argv[2] if len(sys.argv) > 2 else 'kap1'
grid = np.linspace(-0.9, 0.9, ng); t0 = time.time(); ok = bad = 0; badl = []
for R1 in grid:
    for R2 in grid:
        if R2 < R1 - 1e-12: continue
        res = dp(R1, R2, mode)
        if res.status == 0: ok += 1
        else: bad += 1; badl.append((round(R1, 3), round(R2, 3)))
print(f"mode {mode}: shapes {ok + bad}, decomposed-positive certificate found {ok}, not found {bad}  ({time.time() - t0:.0f}s)")
print("not found:", badl[:60])
