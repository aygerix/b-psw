# ak_cover4.py -- admissibility-enforced searches.
# mode 'full': score = min_{P, zeta} max_{Q != P} -log|f_Q(zeta)|        (> 0: fully covered)
# mode 'wta'/'wtal': same with ANY-sign weights theta in [0, 2pi) (score > 0: 0 in conv{(x_Q,y_Q)} at every member);
#   'wtl'/'wtal' include the linear member 2t+1 as an extra coverer.
# mode 'wt'  : score = min_P [ - max_{theta in [0, pi/2]} min_{Q != P} (cos th x_Q + sin th y_Q) ],
#              x_Q = log|f_Q(w_P)|, y_Q = log|f_Q(v_P)|          (> 0: every member unsafe for EVERY nonnegative weighting)
# usage: python3 ak_cover4.py mode n AMIN AMAX seed maxiter
import sys, numpy as np
from scipy.optimize import differential_evolution, minimize
from ak_core import logmat, special
mode = sys.argv[1]; n = int(sys.argv[2]); AMIN = float(sys.argv[3]); AMAX = float(sys.argv[4]); seed = int(sys.argv[5]); MI = int(sys.argv[6])
TH = np.linspace(0, np.pi/2, 181) if mode in ('wt', 'wtl') else np.linspace(0, 2*np.pi, 721)[:-1]
def ac(x):
    a = x[:n]; c = a/2 + (np.sqrt(8*a) - a/2)*x[n:]; return a, c
def score(x):
    if np.any(x[:n] < AMIN) or np.any(x[:n] > AMAX) or np.any(np.abs(x[n:]) >= 0.9999): return -10.0
    a, c = ac(x); Fw, Fv = logmat(a, c)
    X = np.log(np.abs(Fw)); Y = np.log(np.abs(Fv))        # X[Q, P] = log|f_Q(w_P)|
    if mode == 'full':
        Lw = -X.copy(); Lv = -Y.copy(); np.fill_diagonal(Lw, -np.inf); np.fill_diagonal(Lv, -np.inf)
        s = min(Lw.max(0).min(), Lv.max(0).min())
    else:
        s = np.inf
        if mode in ('wtl', 'wtal'):   # linear member 2t+1 included as an extra coverer: log|2w+1|, log|2v+1|
            w, v = special(a, c); xl = np.log(np.abs(2*w + 1)); yl = np.log(np.abs(2*v + 1))
        for P in range(n):
            m = np.arange(n) != P
            xs = X[m, P]; ys = Y[m, P]
            if mode in ('wtl', 'wtal'): xs = np.append(xs, xl[P]); ys = np.append(ys, yl[P])
            comb = np.cos(TH)[:, None]*xs[None, :] + np.sin(TH)[:, None]*ys[None, :]
            s = min(s, -comb.min(1).max())
    return s if np.isfinite(s) else -10.0
bounds = [(AMIN, AMAX)]*n + [(-0.9998, 0.9998)]*n
r = differential_evolution(lambda x: -score(x), bounds, seed=seed, maxiter=MI, popsize=30, tol=1e-12, polish=False, mutation=(0.5, 1.0), recombination=0.9)
x = r.x
p = minimize(lambda x: -score(x), x, method='Nelder-Mead', options={'maxiter': 20000, 'xatol': 1e-12, 'fatol': 1e-14})
if -p.fun > -r.fun: x = p.x
a, c = ac(x)
print("mode=%s n=%d a in [%g,%g] seed=%d: best score %.6g  a=%s c=%s" % (mode, n, AMIN, AMAX, seed, score(x), np.array2string(a, precision=7, max_line_width=300), np.array2string(c, precision=7, max_line_width=300)), flush=True)
