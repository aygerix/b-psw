# audC_rww.py -- adversarial test of the residual claims (i = smaller-a endpoint, j = larger-a endpoint):
#  R-WW: i, j in W, E^v_i(Q_j) < 0  =>  |f_{Q_i}(w_j)| <= 1   (INTjw; O in E^w_j holds since j in W)
#  R-BW: i in B, j in W             =>  |f_{Q_i}(w_j)| <= 1
#  R-BB: i, j in B                  =>  max over the 4 INT tests >= 0 (Q_other in the closed ellipse; O inside automatically)
# Maximises the violation with Nelder-Mead from many starts (penalty outside the case).
import sys, numpy as np
from scipy.optimize import minimize
rng = np.random.default_rng(int(sys.argv[2]))
def fv(a, c, t): return 2*t + 1 + t*(t+1)*(a*t + c)
def up(A, B, C):
    r = np.roots([A, B, C]); return r[np.argmax(r.imag)]
def inD(a, c): return 0 < a < 32 and c*c < 8*a and (a-c)**2 < 8*a
def region(a, c): return "NWVB"[(a - c > 4) + 2*(c > 4)]
def pt(x, y):
    a = 32/(1+np.exp(-x)); lo, hi = a - np.sqrt(8*a), np.sqrt(8*a); return a, (lo+hi)/2 + (hi-lo)/2*np.tanh(y)
def wv(Q): a, c = Q; return up(a, a + c, c + 2), up(a, c, 2)
case = sys.argv[1]
def viol(x):
    Q1, Q2 = pt(x[0], x[1]), pt(x[2], x[3])
    Qi, Qj = (Q1, Q2) if Q1[0] <= Q2[0] else (Q2, Q1)
    if abs(Qi[0]-Qj[0]) + abs(Qi[1]-Qj[1]) < 1e-12: return -10
    ri, rj = region(*Qi), region(*Qj)
    wi, vi = wv(Qi); wj, vj = wv(Qj)
    if case == "WW":
        if ri != "W" or rj != "W": return -10
        if (-fv(Qj[0], Qj[1], vi)).real - 1 >= 0: return -10
        return np.log(abs(fv(Qi[0], Qi[1], wj)))
    if case == "BW":
        if ri != "B" or rj != "W": return -10
        return np.log(abs(fv(Qi[0], Qi[1], wj)))
    if case == "BB":
        if ri != "B" or rj != "B": return -10
        return min(np.log(abs(fv(Qi[0], Qi[1], wj))), np.log(abs(fv(Qi[0], Qi[1], vj))), np.log(abs(fv(Qj[0], Qj[1], wi))), np.log(abs(fv(Qj[0], Qj[1], vi))))
best = (-np.inf, None); tries = 0
while tries < int(sys.argv[3]):
    x0 = rng.normal(0, 2.5, 4)
    if viol(x0) <= -10: continue
    tries += 1
    r = minimize(lambda x: -viol(x)/(1e-300 + 1), x0, method='Nelder-Mead', options={'xatol': 1e-12, 'fatol': 1e-15, 'maxiter': 6000})
    if -r.fun > best[0]: best = (-r.fun, r.x)
x = best[1]; Q1, Q2 = pt(x[0], x[1]), pt(x[2], x[3])
print("%s: starts %d, max violation (>0 refutes) = %.6g at %s, %s; regions %s %s" % (case, tries, best[0], Q1, Q2, region(*Q1), region(*Q2)))
