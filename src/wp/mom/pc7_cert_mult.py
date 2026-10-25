"""Certificates with w-dependent multipliers: find coefficients beta_{i,g} with sum_i sum_g beta_{i,g} g(w) D_i(w) > 0 on samples,
g in a multiplier basis (monomials in power sums).  Usage: python3 pc7_cert_mult.py nconf basis r1,r2 [r1,r2 ...]
basis: 'const' | 'p1' (1, p1) | 'p12' (1, p1, p2) | 'deg2' (1, p1, p2, p1^2) | 'rat1' (1, p1/(1+p1)) | 'ratk' (1/(1+p1)^k-type)"""
import sys, numpy as np
from scipy.optimize import linprog
from pc7_cert_dlib import load, Evec
at = load(); rng = np.random.default_rng(5); nconf = int(sys.argv[1]); basis = sys.argv[2]
CF = [np.exp(rng.uniform(-14, 14, 7)) for _ in range(nconf)]
for _ in range(nconf // 2):
    m = rng.integers(1, 4); vals = np.exp(rng.uniform(-14, 14, m)); idx = rng.integers(0, m, 7)
    CF.append(vals[idx] * np.exp(rng.normal(0, 1e-3, 7)))
def mults(ws):
    p1, p2 = np.sum(ws), np.sum(ws**2)
    if basis == 'const': return [1.0]
    if basis == 'p1': return [1.0, p1]
    if basis == 'p12': return [1.0, p1, p2]
    if basis == 'deg2': return [1.0, p1, p2, p1 * p1]
    if basis == 'rat1': return [1.0, p1 / (1 + p1)]
    if basis == 'rat3': return [1.0, p1 / (1 + p1), (p1 / (1 + p1))**2, (p1 / (1 + p1))**3]
    if basis == 'scale': return [1.0 / (1 + p1)**k for k in range(0, 10)]
for arg in sys.argv[3:]:
    r1, r2 = map(float, arg.split(','))
    K, A, C = at(r1, r2)
    rows = []
    for c in CF:
        D = Evec(K, A, C, c); g = mults(c); rows.append(np.concatenate([gi * D for gi in g]))
    R = np.array(rows); Rn = R / np.abs(R).sum(axis=1, keepdims=True); d = R.shape[1]
    r = linprog(np.concatenate([np.zeros(d), [-1]]), A_ub=np.hstack([-Rn, np.ones((len(Rn), 1))]), b_ub=np.zeros(len(Rn)),
                bounds=[(-1, 1)] * d + [(None, None)], method='highs')
    print(f"r = ({r1:+.3f},{r2:+.3f}) basis {basis}: margin {r.x[-1]:+.3e}", flush=True)
