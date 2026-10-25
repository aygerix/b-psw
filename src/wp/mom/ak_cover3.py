# ak_cover3.py -- differential evolution for fully covered configurations with AMIN <= a <= AMAX.
# score = min over (member, upper special point) of max_{other} -log|f_other(point)|  (> 0  <=> fully covered)
# Members parametrised by (a, s) with a in [AMIN, AMAX], c = a/2 + (sqrt(8a) - a/2) * s, s in [-0.999, 0.999].
# usage: python3 ak_cover3.py n AMIN AMAX seed maxiter
import sys, numpy as np
from scipy.optimize import differential_evolution, minimize
from ak_core import logmat
n = int(sys.argv[1]); AMIN = float(sys.argv[2]); AMAX = float(sys.argv[3]); seed = int(sys.argv[4]); MI = int(sys.argv[5])
def ac(x):
    a = x[:n]; c = a/2 + (np.sqrt(8*a) - a/2)*x[n:]; return a, c
def score(x):
    a, c = ac(x); Fw, Fv = logmat(a, c)
    Lw = -np.log(np.abs(Fw)); Lv = -np.log(np.abs(Fv))
    np.fill_diagonal(Lw, -np.inf); np.fill_diagonal(Lv, -np.inf)
    s = min(Lw.max(0).min(), Lv.max(0).min())
    return s if np.isfinite(s) else -10.0
bounds = [(AMIN, AMAX)]*n + [(-0.999, 0.999)]*n
r = differential_evolution(lambda x: -score(x), bounds, seed=seed, maxiter=MI, popsize=30, tol=1e-12, polish=False, mutation=(0.5, 1.0), recombination=0.9)
x = r.x
p = minimize(lambda x: -score(x), x, method='Nelder-Mead', options={'maxiter': 20000, 'xatol': 1e-12, 'fatol': 1e-14})
if -p.fun > -r.fun: x = p.x
a, c = ac(x)
print("n=%d a in [%g,%g] seed=%d: best score %.6g  a=%s c=%s" % (n, AMIN, AMAX, seed, score(x), np.array2string(a, precision=6), np.array2string(c, precision=6)), flush=True)
