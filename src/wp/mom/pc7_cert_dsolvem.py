"""Mixed-sign version of pc7_cert_dsolve.py: nneg of the w_j negative (fixed sign pattern). Usage: starts seed nneg
Equations scaled row-wise by the sum of absolute values of their terms (so w -> 0 / infinity is not an artificial attractor).
Usage: python3 pc7_cert_dsolve.py starts seed"""
import sys, numpy as np
from scipy.optimize import least_squares
from pc7_cert_dlib import load, Evec
at = load(); rng = np.random.default_rng(int(sys.argv[2])); k = 7
nneg = int(sys.argv[3]); SGN = np.array([-1.0] * nneg + [1.0] * (k - nneg))   # fixed sign pattern of the w_j
def scaled(z):
    r1, r2 = np.tanh(z[0]), np.tanh(z[1]); ws = SGN * np.exp(z[2:])
    K, A, C = at(r1, r2); E = Evec(K, A, C, ws)
    # scale: |kap| + sum |a| + ... approximated by evaluating with absolute coefficients
    Ea = Evec(np.abs(K), np.abs(A), np.abs(C), np.abs(ws)); return E / Ea
hits = []
for s in range(int(sys.argv[1])):
    z0 = np.concatenate([rng.uniform(-2, 2, 2), rng.uniform(-6, 6, k)])
    try:
        r = least_squares(scaled, z0, method='trf', max_nfev=600, xtol=1e-15, ftol=1e-15)
    except Exception as e:
        continue
    res = np.max(np.abs(scaled(r.x)))
    r1, r2 = np.tanh(r.x[0]), np.tanh(r.x[1])
    hits.append((res, r1, r2, np.sort(SGN * np.exp(r.x[2:]))))
    if res < 1e-6: print(f"  start {s}: res {res:.3e}  r = ({r1:+.6f}, {r2:+.6f})  w = {np.array2string(np.sort(SGN * np.exp(r.x[2:])), precision=4)}", flush=True)
hits.sort(key=lambda h: h[0])
print(f"starts {sys.argv[1]}: best scaled residuals")
for h in hits[:12]:
    print(f"  res {h[0]:.3e}  r = ({h[1]:+.6f}, {h[2]:+.6f})  w = {np.array2string(h[3], precision=4)}")
