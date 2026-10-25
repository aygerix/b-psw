# ak_cover.py -- do FULLY COVERED configurations exist?  (every upper special point w_j, v_j of every member j lies in
# the open sublevel set {|f_i| < 1} of some other member i).  If yes, no single-point Theorem-10 extremal argument can work.
# Objective: maximise  min_{j, zeta in {w_j, v_j}}  max_{i != j}  -log|f_i(zeta)|   (> 0  <=> fully covered).
# usage: python3 ak_cover.py n restarts seed
import sys, numpy as np
from scipy.optimize import minimize
from ak_core import ac_from, special, fval, logmat
n = int(sys.argv[1]); R = int(sys.argv[2]); rng = np.random.default_rng(int(sys.argv[3]))
def score(x):
    a, c = ac_from(x[:n], x[n:])
    Fw, Fv = logmat(a, c)
    Lw = -np.log(np.abs(Fw)); Lv = -np.log(np.abs(Fv))
    np.fill_diagonal(Lw, -np.inf); np.fill_diagonal(Lv, -np.inf)
    return min(Lw.max(0).min(), Lv.max(0).min())
best = (-np.inf, None)
for r in range(R):
    x0 = np.concatenate([rng.normal(0, 2.5, n), rng.normal(0, 1.2, n)])
    res = minimize(lambda x: -score(x), x0, method='Nelder-Mead', options={'maxiter': 6000, 'xatol': 1e-10, 'fatol': 1e-12})
    for rep in range(3):
        res = minimize(lambda x: -score(x), res.x, method='Nelder-Mead', options={'maxiter': 6000, 'xatol': 1e-10, 'fatol': 1e-12})
    sc = -res.fun
    if sc > best[0]:
        best = (sc, res.x)
        a, c = ac_from(res.x[:n], res.x[n:])
        print("r=%d best score %.6g  a=%s c=%s" % (r, sc, np.array2string(a, precision=5), np.array2string(c, precision=5)), flush=True)
print("n=%d: max covering score %.6g (>0 means fully covered configuration exists)" % (n, best[0]))
