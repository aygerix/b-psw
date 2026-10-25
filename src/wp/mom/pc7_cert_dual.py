"""Obstructions to a constant-alpha certificate (D-basis): the dual of the max-margin LP gives configurations w^(k) and
weights theta_k >= 0 with sum theta_k D(w^(k))/|D(w^(k))| ~ 0.  Usage: python3 pc7_cert_dual.py nconf r1,r2"""
import sys, numpy as np
from scipy.optimize import linprog
from pc7_cert_dlib import load, Evec, LAM3
at = load(); rng = np.random.default_rng(5); nconf = int(sys.argv[1])
CF = [np.exp(rng.uniform(-14, 14, 7)) for _ in range(nconf)]
for _ in range(nconf // 2):
    m = rng.integers(1, 4); vals = np.exp(rng.uniform(-14, 14, m)); idx = rng.integers(0, m, 7)
    CF.append(vals[idx] * np.exp(rng.normal(0, 1e-3, 7)))
r1, r2 = map(float, sys.argv[2].split(','))
K, A, C = at(r1, r2)
R = np.array([Evec(K, A, C, c) for c in CF]); Rn = R / np.abs(R).sum(axis=1, keepdims=True)
# feasibility LP: Rn alpha >= 1e-9 * ... use max margin with |alpha|_inf <= 1
Aub = np.hstack([-Rn, np.ones((len(Rn), 1))])
r = linprog(np.concatenate([np.zeros(8), [-1]]), A_ub=Aub, b_ub=np.zeros(len(Aub)), bounds=[(-1, 1)] * 8 + [(None, None)], method='highs')
print("margin", r.x[-1])
du = -r.ineqlin.marginals; idx = np.argsort(-du)[:12]
print("rows with largest dual weight (theta, sorted w, normalised D):")
for i in idx:
    if du[i] <= 1e-12: break
    print(f"  theta {du[i]:.3e}  w = {np.array2string(np.sort(CF[i]), precision=3)}\n      D/|D| = {np.array2string(Rn[i], precision=4)}")
