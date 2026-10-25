"""Targeted search for positive solutions of the D-system near a shape where constant-alpha certificates fail:
starts = the obstructing configurations of the max-margin LP dual (and perturbations); unknowns: r2 and the 7 w_j (r1 fixed),
8 equations scaled by the absolute-value evaluation.  Usage: python3 pc7_cert_tsolve.py nconf r1,r2 nstarts"""
import sys, numpy as np
from scipy.optimize import linprog, least_squares
from pc7_cert_dlib import load, Evec
at = load(); rng = np.random.default_rng(11); nconf = int(sys.argv[1]); r1, r2 = map(float, sys.argv[2].split(',')); ns = int(sys.argv[3])
CF = [np.exp(rng.uniform(-14, 14, 7)) for _ in range(nconf)]
for _ in range(nconf // 2):
    m = rng.integers(1, 4); vals = np.exp(rng.uniform(-14, 14, m)); idx = rng.integers(0, m, 7)
    CF.append(vals[idx] * np.exp(rng.normal(0, 1e-3, 7)))
K, A, C = at(r1, r2)
R = np.array([Evec(K, A, C, c) for c in CF]); Rn = R / np.abs(R).sum(axis=1, keepdims=True)
r = linprog(np.concatenate([np.zeros(8), [-1]]), A_ub=np.hstack([-Rn, np.ones((len(Rn), 1))]), b_ub=np.zeros(len(Rn)), bounds=[(-1, 1)] * 8 + [(None, None)], method='highs')
du = -r.ineqlin.marginals; order = np.argsort(-du)
starts = [CF[i] for i in order[:12] if du[i] > 1e-12]
def scaled(z, R1):
    R2 = np.tanh(z[0]); ws = np.exp(z[1:]); K, A, C = at(R1, R2)
    return Evec(K, A, C, ws) / Evec(np.abs(K), np.abs(A), np.abs(C), ws)
best = []
for s in range(ns):
    w0 = starts[s % len(starts)] * np.exp(rng.normal(0, 0.3 if s >= len(starts) else 0, 7))
    z0 = np.concatenate([[np.arctanh(r2)], np.log(w0)])
    try: res = least_squares(scaled, z0, args=(r1,), method='lm', max_nfev=2000, xtol=1e-15, ftol=1e-15)
    except Exception as e: continue
    f = np.max(np.abs(scaled(res.x, r1))); best.append((f, np.tanh(res.x[0]), np.sort(np.exp(res.x[1:]))))
    print(f"  start {s}: scaled residual {f:.3e}  r2 = {np.tanh(res.x[0]):+.6f}  w = {np.array2string(np.sort(np.exp(res.x[1:])), precision=4)}", flush=True)
