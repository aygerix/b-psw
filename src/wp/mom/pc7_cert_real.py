"""Reality-restricted LP: actual (2,2) pencils need K-_j = c_j H + L and K+_j = c_j G + L' without real roots, i.e. every
w_j = 1/c_j lies in I = {w > 0: (w lam - e1)^2 < 4(e2 + w mu), w^2 lam^2 - 4 w (mu + lam e1) + 4 < 0}.
Constant-alpha max margin with configurations sampled in I^7.  Usage: python3 pc7_cert_real.py nconf r1,r2 ..."""
import sys, numpy as np
from scipy.optimize import linprog
from pc7_cert_dlib import load, Evec
at = load(); rng = np.random.default_rng(5); nconf = int(sys.argv[1])
def interval(r1, r2):
    e1, e2 = r1 + r2, r1 * r2; Dl = (1 - r1**2) * (1 - r2**2); lam, mu = 2 * e1 / Dl, 2 * (1 - e1**2 + e2) / Dl
    ws = np.exp(np.linspace(-30, 30, 200001))
    ok = ((ws * lam - e1)**2 < 4 * (e2 + ws * mu)) & (ws**2 * lam**2 - 4 * ws * (mu + lam * e1) + 4 < 0)
    if not ok.any(): return None
    return ws[ok].min(), ws[ok].max(), ok.mean()
def margin(R):
    Rn = R / np.abs(R).sum(axis=1, keepdims=True); d = R.shape[1]
    r = linprog(np.concatenate([np.zeros(d), [-1]]), A_ub=np.hstack([-Rn, np.ones((len(Rn), 1))]), b_ub=np.zeros(len(Rn)), bounds=[(-1, 1)] * d + [(None, None)], method='highs')
    return r.x[-1] if r.status == 0 else np.nan
for arg in sys.argv[2:]:
    r1, r2 = map(float, arg.split(',')); I = interval(r1, r2)
    if I is None: print(f"r = ({r1:+.3f},{r2:+.3f}): reality interval EMPTY -> no real pencil at this shape"); continue
    lo, hi, _ = I; L0, L1 = np.log(lo), np.log(hi)
    CF = [np.exp(rng.uniform(L0, L1, 7)) for _ in range(nconf)]
    for _ in range(nconf // 2):
        m = rng.integers(1, 4); vals = np.exp(rng.uniform(L0, L1, m)); idx = rng.integers(0, m, 7)
        CF.append(np.clip(vals[idx] * np.exp(rng.normal(0, 1e-3, 7)), lo, hi))
    K, A, C = at(r1, r2); D = np.array([Evec(K, A, C, c) for c in CF])
    print(f"r = ({r1:+.3f},{r2:+.3f}): reality interval w in [{lo:.4g}, {hi:.4g}]  constant-alpha margin {margin(D):+.3e}", flush=True)
