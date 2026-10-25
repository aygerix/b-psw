# audC_rob.py -- extension of C2 to collinear cubics + 2t+1: ROBUST endpoint tests (independent of intermediate members).
# Endpoints Q_i, Q_j in D (segment [Q_i,Q_j] in D by convexity), linear member O = (0,0).
# For endpoint i (other endpoint j) and z in {w_i, v_i}: factor map F_z(Q) = f_Q(w_i) (z = w) or -f_Q(v_i) (z = v);
# cubic factors on the segment are 1 + s Z, s in (0,1], Z = F_z(Q_j) - 1;  O-factor F_z(O).
#   EXT(i,z): Re Z >= 0 and |F_z(O)| >= 1      (every cubic on the open segment has |factor| > 1)
#   INT(i,z): |F_z(Q_j)| <= 1 and |F_z(O)| <= 1 (every cubic on the segment has |factor| <= 1, interior ones < 1)
# score of a test = min of its normalised conditions (>= 0 iff it passes); C = max over the 8 tests; C < 0 = not robustly covered.
import sys, numpy as np
from scipy.optimize import minimize
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 3)
def fv(a, c, t): return 2*t + 1 + t*(t+1)*(a*t + c)
def up(A, B, C):
    r = np.roots([A, B, C]); return r[np.argmax(r.imag)]
def inD(a, c): return a > 0 and c*c < 8*a and (a-c)**2 < 8*a
def pt(x, y):
    a = 32/(1 + np.exp(-x)); lo = a - np.sqrt(8*a); hi = np.sqrt(8*a)
    return a, (lo + hi)/2 + (hi - lo)/2*np.tanh(y)
def tests(Qi, Qj):
    out = {}
    for tag, (A, B) in (("i", (Qi, Qj)), ("j", (Qj, Qi))):
        a, c = A; w = up(a, a + c, c + 2); v = up(a, c, 2)
        for zn, F in (("w", lambda Q: fv(Q[0], Q[1], w)), ("v", lambda Q: -fv(Q[0], Q[1], v))):
            Z = F(B) - 1; FO = F((0.0, 0.0))
            out["EXT" + tag + zn] = min(Z.real/abs(Z), np.log(abs(FO)))
            out["INT" + tag + zn] = min(-np.log(abs(1 + Z)), -np.log(abs(FO)))
    return out
def Cscore(x):
    Qi, Qj = pt(x[0], x[1]), pt(x[2], x[3])
    if not (inD(*Qi) and inD(*Qj)) or np.hypot(Qi[0]-Qj[0], Qi[1]-Qj[1]) < 1e-9: return 10.0, None
    t = tests(Qi, Qj); return max(t.values()), t
mode = sys.argv[1]
if mode == "rand":
    T = int(sys.argv[3]); unc = 0; worst = (np.inf, None); passcount = {}
    for _ in range(T):
        x = rng.normal(0, 2.5, 4); C, t = Cscore(x)
        if t is None: continue
        if C < 0: unc += 1
        if C < worst[0]: worst = (C, x)
        for kk, vv in t.items():
            if vv >= 0: passcount[kk] = passcount.get(kk, 0) + 1
    print("rand pairs %d: not robustly covered %d; min C = %.4g at %s" % (T, unc, worst[0], [pt(worst[1][0], worst[1][1]), pt(worst[1][2], worst[1][3])]))
    print("pass counts:", dict(sorted(passcount.items())))
else:
    S = int(sys.argv[3]); best = (np.inf, None)
    for it in range(S):
        x0 = rng.normal(0, 2.5, 4)
        if Cscore(x0)[1] is None: continue
        r = minimize(lambda x: Cscore(x)[0], x0, method='Nelder-Mead', options={'xatol': 1e-11, 'fatol': 1e-13, 'maxiter': 8000})
        if r.fun < best[0]: best = (r.fun, r.x)
    x = best[1]; C, t = Cscore(x)
    print("opt starts %d: min C = %.6g" % (S, C)); print("Qi, Qj =", pt(x[0], x[1]), pt(x[2], x[3])); print("tests:", t)
