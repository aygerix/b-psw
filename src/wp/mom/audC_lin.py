# audC_lin.py -- audit of line ak, extension of C2: cubic members with collinear (a,c) points PLUS 2t+1.
# New code (does not import ak_*). Members f = 2t+1 + t(t+1)(a t + c) with (a,c) in
# D = {a > 0, c^2 < 8a, (a-c)^2 < 8a}; linear member 2t+1 = point (0,0).
# For a configuration, at each cubic member i and each upper special point z in {w_i, v_i}:
#   L_l = log|f_l(z)| over l != i (including the linear member).
#   s(i,z) = min(max_l L_l, -min_l L_l)  (> 0 iff the moduli straddle 1, i.e. the single-point test fails at z).
# (cubic logs are divided by |Q_l - Q_i| so that coalescing members do not drive J to 0; sign pattern unchanged)
# J = min_{i,z} s(i,z) > 0 means: NO single special point gives a contradiction ("uncovered").
# Mode 'rand': random collinear configs; mode 'opt': maximise J (Nelder-Mead from random starts).
import sys, numpy as np
from scipy.optimize import minimize
rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
def fv(a, c, t): return 2*t + 1 + t*(t+1)*(a*t + c)
def up(A, B, C):
    r = np.roots([A, B, C]); return r[np.argmax(r.imag)]
def inD(a, c): return a > 0 and c*c < 8*a and (a-c)**2 < 8*a
def pt(x, y):
    # map R^2 -> D: a = exp(x), c = a/2 + (sqrt(8a) - a/2)*tanh(y) if a < 32 ... use the exact interval of c
    a = np.exp(x)
    # c-interval of D at fixed a: max(-sqrt(8a), a - sqrt(8a)) < c < min(sqrt(8a), a + sqrt(8a))
    lo = max(-np.sqrt(8*a), a - np.sqrt(8*a)); hi = min(np.sqrt(8*a), a + np.sqrt(8*a))
    if hi <= lo: return None
    return a, (lo + hi)/2 + (hi - lo)/2*np.tanh(y)
GEN = False
ONESIDED = "onesided" in sys.argv   # calibration: s = -min L > 0 iff the one-sided (all >= 1) test fails (C6 'fully covered')
def config(x, m):
    if GEN:
        P = [pt(x[2*j], x[2*j+1]) for j in range(m)]
        return None if any(q is None for q in P) else P
    # x = [xA, yA, xB, yB, tau_1..tau_{m-2}]  -> m points on the segment [QA, QB]
    PA = pt(x[0], x[1]); PB = pt(x[2], x[3])
    if PA is None or PB is None: return None
    taus = [0.0, 1.0] + [1/(1+np.exp(-z)) for z in x[4:4+m-2]]
    return [(PA[0] + s*(PB[0]-PA[0]), PA[1] + s*(PB[1]-PA[1])) for s in taus]
def J(cfg, lin=True):
    m = len(cfg); best = np.inf; detail = None
    for i, (a, c) in enumerate(cfg):
        if not inD(a, c): return -np.inf, None
        w = up(a, a + c, c + 2); v = up(a, c, 2)
        for name, z in (("w", w), ("v", v)):
            L = [np.log(abs(fv(al, cl, z)))/np.hypot(al - a, cl - c) for l, (al, cl) in enumerate(cfg) if l != i]
            if lin: L.append(np.log(abs(2*z + 1)))
            s = -min(L) if ONESIDED else min(max(L), -min(L))
            if s < best: best, detail = s, (i, name, L)
    return best, detail
def negJ(x, m, lin):
    cfg = config(x, m)
    if cfg is None: return 10.0
    j, _ = J(cfg, lin)
    return 10.0 if not np.isfinite(j) else -j
mode = sys.argv[1]; m = int(sys.argv[2]);
if mode.startswith("gen"): GEN = True; mode = mode[3:]
def x0gen(m):
    if GEN: return np.concatenate([[rng.normal(0, 2.5), rng.normal(0, 1.5)] for _ in range(m)])
    return np.concatenate([rng.normal(0, 2.5, 1), rng.normal(0, 1.5, 1), rng.normal(0, 2.5, 1), rng.normal(0, 1.5, 1), rng.normal(0, 2, max(m-2, 0))])
lin = (len(sys.argv) <= 4 or sys.argv[4] != "nolin")
if mode == "rand":
    T = int(sys.argv[5]) if len(sys.argv) > 5 else 20000
    cnt = 0; top = -np.inf; topx = None
    for _ in range(T):
        x = x0gen(m)
        cfg = config(x, m)
        if cfg is None: continue
        j, d = J(cfg, lin)
        if j > 0: cnt += 1
        if j > top: top, topx = j, x
    print("rand m=%d lin=%s: uncovered %d/%d, max J = %.4g" % (m, lin, cnt, T, top)); print("argmax x =", list(topx)); print("cfg =", config(topx, m))
else:
    S = int(sys.argv[5]) if len(sys.argv) > 5 else 40
    best = (np.inf, None)
    for it in range(S):
        x0 = x0gen(m)
        if negJ(x0, m, lin) >= 10: continue
        r = minimize(negJ, x0, args=(m, lin), method='Nelder-Mead', options={'xatol': 1e-10, 'fatol': 1e-13, 'maxiter': 6000})
        if r.fun < best[0]: best = (r.fun, r.x)
    x = best[1]; cfg = config(x, m); j, d = J(cfg, lin)
    print("opt m=%d lin=%s starts=%d: max J = %.6g" % (m, lin, S, j)); print("x =", list(x)); print("cfg =", cfg); print("worst (i, point, logs) =", d)
