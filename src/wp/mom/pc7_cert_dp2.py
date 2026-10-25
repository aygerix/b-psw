"""Relaxation ladder at fixed shapes (all sufficient for sum alpha_i E_i > 0 on (0, inf)^7):
  dp   : kap = 1, a >= 0, b >= 0, lam3 >= 0
  pair : kap = 1, pi(w, v) = kap/21 + (a(w) + a(v))/6 + b(w, v) >= 0, lam3 >= 0
  trip : kap = 1, tau(u,v,w) = kap/35 + (a(u)+a(v)+a(w))/15 + (b(u,v)+b(u,w)+b(v,w))/5 + lam3 u v w >= 0
  full : kap free, E_alpha(w) >= 1e-6 on random 7-member samples (the original LP, for comparison)
Usage: python3 pc7_cert_dp2.py ngrid modes"""
import sys, time, itertools, numpy as np
from scipy.optimize import linprog
from pc7_cert_lib import load, at_shape, LAM3
D = load(); rng = np.random.default_rng(1)
WG = np.exp(np.linspace(np.log(1e-6), np.log(1e6), 241))
PW = np.exp(np.linspace(np.log(1e-5), np.log(1e5), 46))
WV = np.array([(a, b) for a in PW for b in PW if b >= a])
P3 = np.exp(np.linspace(np.log(1e-4), np.log(1e4), 19))
UVW = np.array([(a, b, c) for a in P3 for b in P3 for c in P3 if a <= b <= c])
def aval(Ac, x): return np.array([np.polyval(Ac[i][::-1], x) for i in range(8)])
def bval(Bc, a, b):
    M = np.outer(a ** np.arange(10), b ** np.arange(10)); return np.array([np.sum(Bc[i] * M) for i in range(8)])
def nrm(R): return R / np.maximum(np.abs(R).max(axis=1, keepdims=True), 1e-300)
def solve(Aub, Aeq=None, beq=None):
    return linprog(np.zeros(8), A_ub=-Aub, b_ub=np.zeros(len(Aub)), A_eq=Aeq, b_eq=beq, bounds=[(-1e4, 1e4)] * 8, method='highs').status == 0
def test(R1, R2, mode):
    kap, Ac, Bc = at_shape(D, R1, R2)
    if mode == 'dp':
        R = np.vstack([nrm(np.array([aval(Ac, x) for x in WG])), nrm(np.array([bval(Bc, a, b) for a, b in WV])), LAM3[None, :]])
        return solve(R, kap[None, :], [1.0])
    if mode == 'pair':
        R = np.vstack([nrm(np.array([kap / 21 + (aval(Ac, a) + aval(Ac, b)) / 6 + bval(Bc, a, b) for a, b in WV])), LAM3[None, :]])
        return solve(R, kap[None, :], [1.0])
    if mode == 'trip':
        R = nrm(np.array([kap / 35 + (aval(Ac, a) + aval(Ac, b) + aval(Ac, c)) / 15 + (bval(Bc, a, b) + bval(Bc, a, c) + bval(Bc, b, c)) / 5 + LAM3 * a * b * c for a, b, c in UVW]))
        return solve(R, kap[None, :], [1.0])
    if mode == 'full':
        rows = []
        for _ in range(800):
            ws = np.exp(rng.uniform(-8, 6, 7)); E = kap.copy()
            for x in ws: E = E + aval(Ac, x)
            for x, y in itertools.combinations(ws, 2): E = E + bval(Bc, x, y)
            E = E + LAM3 * sum(np.prod(c) for c in itertools.combinations(ws, 3)); rows.append(E)
        R = nrm(np.array(rows))
        return linprog(np.zeros(8), A_ub=-R, b_ub=-1e-6 * np.ones(len(R)), bounds=[(-1e4, 1e4)] * 8, method='highs').status == 0
ng = int(sys.argv[1]); modes = sys.argv[2].split(',')
grid = np.linspace(-0.9, 0.9, ng)
shapes = [(a, b) for a in grid for b in grid if b >= a - 1e-12]
for mode in modes:
    t0 = time.time(); res = [(s, test(*s, mode)) for s in shapes]
    bad = [(round(s[0], 3), round(s[1], 3)) for s, ok in res if not ok]
    print(f"mode {mode}: shapes {len(res)}, certificate found {len(res) - len(bad)}, not found {len(bad)} ({time.time() - t0:.0f}s)  not found: {bad}", flush=True)
