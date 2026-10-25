# ak_dense.py -- large-n test of WS' (any-sign weighted single-member safety) on dense convex configurations:
# n members on an ellipse in the (V, ut) plane (V = 1/a, ut = c/a - 1/2). For each member P: margin_P =
# max_theta min_{Q != P} (cos th log|f_Q(w_P)| + sin th log|f_Q(v_P)|). Violation of WS' <=> margin_P <= 0 for ALL P.
# Prints the configurations with the largest  min_P(-margin_P)  (> 0 = violation).
import sys, numpy as np
from ak_core import logmat
rng = np.random.default_rng(int(sys.argv[1])); NC = int(sys.argv[2])
TH = np.linspace(0, 2*np.pi, 1441)[:-1]; C, S = np.cos(TH), np.sin(TH)
def admissible(V, ut): return (V > 1/32) & ((np.abs(ut) + 0.5)**2 < 8*V)
best = []
for it in range(NC):
    V0 = np.exp(rng.uniform(np.log(1/32), np.log(4.0))); lim = np.sqrt(8*V0) - 0.5
    u0 = rng.uniform(-1, 1)*lim
    if not admissible(V0, u0): continue
    n = rng.choice([24, 36, 48]); r = np.exp(rng.uniform(np.log(1e-4), np.log(0.3)))*min(V0, 1)
    ecc = np.exp(rng.uniform(-2.5, 2.5)); rot = rng.uniform(0, np.pi)
    ph = np.linspace(0, 2*np.pi, n, endpoint=False) + rng.uniform(0, 1)
    x, y = r*np.cos(ph)*ecc, r*np.sin(ph)/ecc
    V = V0 + x*np.cos(rot) - y*np.sin(rot); U = u0 + x*np.sin(rot) + y*np.cos(rot)
    if not np.all(admissible(V, U)): continue
    a = 1/V; c = a*(U + 0.5)
    Fw, Fv = logmat(a, c); X = np.log(np.abs(Fw)); Y = np.log(np.abs(Fv))
    margins = np.empty(n)
    for P in range(n):
        m = np.arange(n) != P
        comb = C[:, None]*X[m, P][None, :] + S[:, None]*Y[m, P][None, :]
        margins[P] = comb.min(1).max()
    sc = -margins.max()
    best.append((sc, V0, u0, r, ecc, rot, n, int(np.sum(margins <= 0))))
best.sort(reverse=True)
print("tested %d dense configurations; top scores (score, V0, ut0, r, ecc, rot, n, #unsafe members):" % len(best))
for b in best[:8]: print("  ", tuple(float('%.4g' % z) for z in b))
print("violations (all members unsafe):", sum(1 for b in best if b[0] > 0))
