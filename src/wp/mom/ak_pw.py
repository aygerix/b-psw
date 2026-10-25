# ak_pw.py -- test the pairwise conjecture (PW):  Q(i,j) := log|f_i(w_j) f_i(v_j)|  has the sign of a_i - a_j.
# (w_j, v_j: upper "-" and "+" special points of member j.)  Also checks the first-order identities
#   Re p(w) = Re p(v),  Re[w p(w) - v p(v)] = 3/a   (p = t(t+1)).
import sys, numpy as np
from scipy.optimize import minimize
from ak_core import ac_from, special, fval
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
def Q(ai, ci, aj, cj):
    w, v = special(np.atleast_1d(aj), np.atleast_1d(cj))
    return np.log(np.abs(fval(ai, ci, w)*fval(ai, ci, v)))
# identities
al = rng.normal(0, 3, 10000); ga = rng.normal(0, 1.5, 10000); a, c = ac_from(al, ga)
w, v = special(a, c); p = lambda t: t*(t+1)
print("max |Re p(w) - Re p(v)| =", np.max(np.abs(p(w).real - p(v).real)/(1+np.abs(p(w)))))
print("max |a Re[w p(w) - v p(v)] - 3| =", np.max(np.abs(a*(w*p(w) - v*p(v)).real - 3)))
# random pairs
M = 400000
ai, ci = ac_from(rng.normal(0, 3, M), rng.normal(0, 1.5, M)); aj, cj = ac_from(rng.normal(0, 3, M), rng.normal(0, 1.5, M))
q = Q(ai, ci, aj, cj); s = np.sign(ai - aj)
bad_up = np.sum((s > 0) & (q <= 0)); bad_dn = np.sum((s < 0) & (q >= 0))
print("random pairs: violations a_i>a_j & Q<=0: %d ;  a_i<a_j & Q>=0: %d  (of %d)" % (bad_up, bad_dn, M))
# adversarial: minimise Q subject to a_i > a_j (penalised), and maximise Q subject to a_i < a_j
def obj_up(x):
    ai, ci = ac_from(x[0], x[1]); aj, cj = ac_from(x[2], x[3])
    if ai <= aj: return 10.0 + (aj - ai)
    return float(Q(ai, ci, aj, cj)[0])/ (np.log(ai/aj))     # normalised by log(a_i/a_j)
def obj_dn(x):
    ai, ci = ac_from(x[0], x[1]); aj, cj = ac_from(x[2], x[3])
    if ai >= aj: return 10.0 + (ai - aj)
    return -float(Q(ai, ci, aj, cj)[0])/ (np.log(aj/ai))
for name, ob in (("up", obj_up), ("dn", obj_dn)):
    best = (1e9, None)
    for it in range(60):
        x0 = np.concatenate([rng.normal(0, 3, 1), rng.normal(0, 1.5, 1), rng.normal(0, 3, 1), rng.normal(0, 1.5, 1)])
        r = minimize(ob, x0, method='Nelder-Mead', options={'xatol': 1e-12, 'fatol': 1e-14, 'maxiter': 20000})
        if r.fun < best[0]: best = (r.fun, r.x)
    x = best[1]; ai, ci = ac_from(x[0], x[1]); aj, cj = ac_from(x[2], x[3])
    print("adversarial %s: min normalised value %.6g at a_i=%.6g c_i=%.6g a_j=%.6g c_j=%.6g  Q=%.3g" % (name, best[0], ai, ci, aj, cj, float(Q(ai, ci, aj, cj)[0])))
