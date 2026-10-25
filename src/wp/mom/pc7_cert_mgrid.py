"""Grid survey over the shape triangle (e1 >= 0 suffices by t -> -t): max margin of constant-alpha and {1, p1}-multiplier
certificates, D-basis, broad sampling.  Usage: python3 pc7_cert_mgrid.py nconf step"""
import sys, numpy as np
from scipy.optimize import linprog
from pc7_cert_dlib import load, Evec
at = load(); rng = np.random.default_rng(5); nconf = int(sys.argv[1]); step = float(sys.argv[2])
CF = [np.exp(rng.uniform(-14, 14, 7)) for _ in range(nconf)]
for _ in range(nconf // 2):
    m = rng.integers(1, 4); vals = np.exp(rng.uniform(-14, 14, m)); idx = rng.integers(0, m, 7)
    CF.append(vals[idx] * np.exp(rng.normal(0, 1e-3, 7)))
P1 = np.array([np.sum(c) for c in CF])
def margin(R):
    Rn = R / np.abs(R).sum(axis=1, keepdims=True); d = R.shape[1]
    r = linprog(np.concatenate([np.zeros(d), [-1]]), A_ub=np.hstack([-Rn, np.ones((len(Rn), 1))]), b_ub=np.zeros(len(Rn)), bounds=[(-1, 1)] * d + [(None, None)], method='highs')
    return r.x[-1] if r.status == 0 else np.nan
vals = np.round(np.arange(-0.95, 0.951, step), 4)
res = []
for a in vals:
    for b in vals:
        if b < a or a + b < -1e-9: continue
        K, A, C = at(a, b); D = np.array([Evec(K, A, C, c) for c in CF])
        m0 = margin(D); m1 = margin(np.hstack([D, P1[:, None] * D]))
        res.append((a, b, m0, m1)); print(f"r = ({a:+.2f},{b:+.2f})  const {m0:+.2e}   1,p1 {m1:+.2e}", flush=True)
res = np.array(res)
print("SUMMARY: shapes", len(res), " const-alpha margin <= 1e-9 at", int(np.sum(res[:, 2] <= 1e-9)), " {1,p1} margin <= 1e-9 at", int(np.sum(res[:, 3] <= 1e-9)))
print("{1,p1} failures:", [(r[0], r[1]) for r in res if r[3] <= 1e-9])
