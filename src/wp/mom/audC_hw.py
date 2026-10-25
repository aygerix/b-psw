# audC_hw.py -- test hypothesis H-W: for X in W u B (O inside E^w_X, i.e. a_X - c_X > 4) and Y in D with
# a_Y <= a_X and c_Y >= c_X:  |f_Y(w_X)| <= 1  (Y in the closed ellipse E^w_X).  Random + adversarial.
# Also: BB pairs with 0 < c_j - c_i < a_j - a_i: which tests pass.
import sys, numpy as np, collections
from scipy.optimize import minimize
rng = np.random.default_rng(int(sys.argv[1]))
YWB = len(sys.argv) > 4 and sys.argv[4] == "YWB"   # H-W': also require a_Y - c_Y > 4
def fv(a, c, t): return 2*t + 1 + t*(t+1)*(a*t + c)
def up(A, B, C):
    r = np.roots([A, B, C]); return r[np.argmax(r.imag)]
def inD(a, c): return 0 < a < 32 and c*c < 8*a and (a-c)**2 < 8*a
def region(a, c): return "NWVB"[(a - c > 4) + 2*(c > 4)]
def sampD():
    while True:
        a = 32*rng.random(); lo, hi = a - np.sqrt(8*a), np.sqrt(8*a); c = lo + (hi - lo)*rng.random()
        if inD(a, c): return (a, c)
def sampWB():
    while True:
        X = sampD()
        if X[0] - X[1] > 4: return X
viol = 0; n = 0; worst = -np.inf; wx = None
for _ in range(int(sys.argv[2])):
    X = sampWB(); Y = sampD()
    if rng.random() < 0.5:  # bias Y toward the quadrant
        Y = (X[0] - (X[0])*rng.random()**2, X[1] + 6*rng.random()**2)
        if not inD(*Y): continue
    if not (Y[0] <= X[0] and Y[1] >= X[1]): continue
    if YWB and not (Y[0] - Y[1] > 4): continue
    n += 1; w = up(X[0], X[0] + X[1], X[1] + 2); val = np.log(abs(fv(Y[0], Y[1], w)))
    if val > worst: worst, wx = val, (X, Y)
    if val > 1e-12: viol += 1
print("H-W random: n=%d violations=%d max log|f_Y(w_X)| = %.4g at %s" % (n, viol, worst, wx))
# adversarial: maximise log|f_Y(w_X)| over the constraint set (penalised)
def obj(x):
    a = 32/(1+np.exp(-x[0])); lo, hi = a - np.sqrt(8*a), np.sqrt(8*a); c = (lo+hi)/2 + (hi-lo)/2*np.tanh(x[1])
    if a - c <= 4: return 10 + (4 - (a - c))
    s1 = 1/(1+np.exp(-x[2])); aY = a*s1; loY, hiY = aY - np.sqrt(8*aY), np.sqrt(8*aY)
    cY = (loY+hiY)/2 + (hiY-loY)/2*np.tanh(x[3])
    if cY < c: return 10 + (c - cY)
    if YWB and aY - cY <= 4: return 10 + (4 - (aY - cY))
    w = up(a, a + c, c + 2); return -np.log(abs(fv(aY, cY, w)))
best = (np.inf, None)
for it in range(int(sys.argv[3])):
    x0 = rng.normal(0, 2, 4)
    if obj(x0) >= 10: continue
    r = minimize(obj, x0, method='Nelder-Mead', options={'xatol': 1e-12, 'fatol': 1e-14, 'maxiter': 5000})
    if r.fun < best[0]: best = (r.fun, r.x)
print("H-W adversarial: max log|f_Y(w_X)| = %.6g" % (-best[0]))
x = best[1]; a = 32/(1+np.exp(-x[0])); lo, hi = a - np.sqrt(8*a), np.sqrt(8*a); c = (lo+hi)/2 + (hi-lo)/2*np.tanh(x[1])
aY = a/(1+np.exp(-x[2])); loY, hiY = aY - np.sqrt(8*aY), np.sqrt(8*aY); cY = (loY+hiY)/2 + (hiY-loY)/2*np.tanh(x[3])
print("   at X =", (a, c), "Y =", (aY, cY))
