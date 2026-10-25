# ak_kmin.py -- local curvature threshold kappa_min(P, phi) for WS' at a vertex P with unit outer normal N = (cos phi, sin phi)
# (coordinates V = 1/a, ut = c/a - 1/2).  A dense convex polygon along Gamma has P unsafe iff kappa_Gamma(P) < kappa_min.
# kappa_min = a [lam Re(Lw(e)^2) + mu Re(Lv(e)^2)] / |N|, N = lam d^w + mu d^v, d^w = (3/2+ut, rho), d^v = (3/2-ut, -rho),
# e = (-sin phi, cos phi), Lw(e) = 2 e_V w + e_ut p(w), Lv(e) = 2 e_V (v+1) + e_ut p(v).
# Also a brute-force check of the formula: second difference of Lambda along e at P (Q = P + s e) vs the formula.
import numpy as np
from ak_core import special, fval
PH = np.linspace(0, 2*np.pi, 721)[:-1]
def kmin(V, ut):
    a = 1/V; c = a*(ut + 0.5); w, v = special(np.array([a]), np.array([c])); w, v = w[0], v[0]
    rho = 2*V + 1/8 - ut**2/2; mw, mv = 1.5 + ut, 1.5 - ut
    lam = (np.cos(PH) + mv*np.sin(PH)/rho)/3; mu = (np.cos(PH) - mw*np.sin(PH)/rho)/3
    eV, eU = -np.sin(PH), np.cos(PH)
    Lw = 2*eV*w + eU*w*(w+1); Lv = 2*eV*(v+1) + eU*v*(v+1)
    return a*(lam*(Lw**2).real + mu*(Lv**2).real), lam, mu, a, c, w, v
# formula check by finite differences at random points
rng = np.random.default_rng(3); errs = []
for it in range(200):
    V = np.exp(rng.uniform(np.log(1/32), np.log(3))); lim = np.sqrt(8*V) - 0.5; ut = rng.uniform(-0.9, 0.9)*lim
    K, lam, mu, a, c, w, v = kmin(V, ut); j = rng.integers(len(PH)); ph = PH[j]; e = np.array([-np.sin(ph), np.cos(ph)])
    def Lam(s):
        V2, u2 = V + s*e[0], ut + s*e[1]; a2 = 1/V2; c2 = a2*(u2 + 0.5)
        return lam[j]*np.log(abs(fval(a2, c2, w))) + mu[j]*np.log(abs(fval(a2, c2, v)))
    h = 1e-4*min(V, lim + 1e-9)
    d2 = (Lam(h) - 2*Lam(0) + Lam(-h))/h**2
    errs.append(abs(-d2/a - K[j]) / (1 + abs(K[j])))     # e^T H e = -a^2[...]  =>  kappa_min = -d2/a (|N| = 1)
print("formula check: max rel err %.3g (median %.3g)" % (max(errs), np.median(errs)))
# scan the domain: min over phi of kappa_min, and the fraction of bad directions
rows = []
for V in np.exp(np.linspace(np.log(1/32 + 1e-4), np.log(20), 120)):
    lim = np.sqrt(8*V) - 0.5
    for ut in np.linspace(-lim, lim, 41)[1:-1]:
        K = kmin(V, ut)[0]
        rows.append((K.min(), np.mean(K > 0), 1/V, ut))
rows = np.array(rows)
print("points scanned:", len(rows), "; points with ALL directions bad (min_phi kappa_min > 0):", int(np.sum(rows[:, 0] > 0)))
i = np.argmax(rows[:, 0]); print("max over P of min_phi kappa_min = %.4g at a=%.4g ut=%.4g (bad fraction %.3f)" % (rows[i, 0], rows[i, 2], rows[i, 3], rows[i, 1]))
print("bad-direction fraction: min %.3f median %.3f max %.3f" % (rows[:, 1].min(), np.median(rows[:, 1]), rows[:, 1].max()))
