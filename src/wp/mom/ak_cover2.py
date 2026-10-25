# ak_cover2.py -- fully covered configurations with a <= AMAX (away from the Chebyshev corner a = 32).
# usage: python3 ak_cover2.py n restarts seed AMAX
import sys, numpy as np
from scipy.optimize import minimize
from ak_core import special, fval, logmat
n = int(sys.argv[1]); R = int(sys.argv[2]); rng = np.random.default_rng(int(sys.argv[3])); AM = float(sys.argv[4])
def ac(x):
    a = AM/(1.0 + np.exp(-x[:n])); c = a/2 + (np.sqrt(8*a) - a/2)*np.tanh(x[n:]); return a, c
def score(x):
    a, c = ac(x); Fw, Fv = logmat(a, c)
    Lw = -np.log(np.abs(Fw)); Lv = -np.log(np.abs(Fv))
    np.fill_diagonal(Lw, -np.inf); np.fill_diagonal(Lv, -np.inf)
    return min(Lw.max(0).min(), Lv.max(0).min())
best = (-np.inf, None)
for r in range(R):
    x0 = np.concatenate([rng.normal(0, 2.5, n), rng.normal(0, 1.2, n)])
    res = minimize(lambda x: -score(x), x0, method='Nelder-Mead', options={'maxiter': 8000, 'xatol': 1e-11, 'fatol': 1e-13})
    for rep in range(4):
        res = minimize(lambda x: -score(x), res.x, method='Nelder-Mead', options={'maxiter': 8000, 'xatol': 1e-11, 'fatol': 1e-13})
    if -res.fun > best[0]:
        best = (-res.fun, res.x); a, c = ac(res.x)
        print("r=%d score %.6g a=%s c=%s" % (r, -res.fun, np.array2string(a, precision=6), np.array2string(c, precision=6)), flush=True)
print("n=%d AMAX=%g: best covering score %.6g" % (n, AM, best[0]))
