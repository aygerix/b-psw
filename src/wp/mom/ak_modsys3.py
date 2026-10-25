# ak_modsys2.py -- modulus-only system with deflation in the chart P = (log a, c/sqrt a) and a positive control.
# usage: python3 ak_modsys2.py k nstarts seed control(0/1) [lin]
#   control = 1: targets tau = R(X*) for a random X* (a solution exists); measures solver recall.
import sys, numpy as np
from scipy.optimize import least_squares
from ak_core import ac_from, special, fval, logmat
k = int(sys.argv[1]); N = int(sys.argv[2]); seed = int(sys.argv[3]); ctl = int(sys.argv[4])
lin = int(sys.argv[5]) if len(sys.argv) > 5 else 0
rng = np.random.default_rng(seed); m = k - lin
def unpack(x): return ac_from(x[:m], x[m:])
def R(x):
    a, c = unpack(x)
    Fw, Fv = logmat(a, c)
    out = np.concatenate([np.log(np.abs(Fw)).sum(0), np.log(np.abs(Fv)).sum(0)])
    if lin:
        w, v = special(a, c)
        out = out + np.concatenate([np.log(np.abs(2*w + 1)), np.log(np.abs(2*v + 1))])
    return out
def chart(x):
    a, c = unpack(x); return np.stack([np.log(a), c/np.sqrt(a)], 1)
tau = np.zeros(2*m)
if ctl:
    xs = np.concatenate([rng.normal(0, 2.5, m), rng.normal(0, 1.2, m)]); tau = R(xs)
iu = np.triu_indices(m, 1)
def Rd(x):
    P = chart(x); D2 = ((P[:, None, :] - P[None, :, :])**2).sum(-1)[iu]
    return (R(x) - tau)*np.prod(1.0 + 0.05/D2)
nsol = 0; best = []
for it in range(N):
    x0 = np.concatenate([rng.normal(0, 2.5, m), rng.normal(0, 1.2, m)])
    try:
        r = least_squares(Rd, x0, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=1500)
    except Exception:
        continue
    x = r.x
    if not np.all(np.isfinite(x)): continue
    res = np.max(np.abs(R(x) - tau)); P = chart(x)
    D = np.sqrt(((P[:, None, :] - P[None, :, :])**2).sum(-1))[iu]; dmin = D.min() if len(D) else 1.0
    a, c = unpack(x)
    best.append((res, dmin))
    print("it=%d res=%.3e dmin=%.3e" % (it, res, dmin), flush=True)
    if res < 1e-10 and dmin > 1e-3:
        nsol += 1
        if nsol <= 5: print("SOL it=%d res=%.2e dmin=%.3g a=%s c=%s" % (it, res, dmin, np.array2string(a, precision=5), np.array2string(c, precision=5)), flush=True)
best.sort()
print("k=%d lin=%d ctl=%d: %d/%d genuine; 5 smallest (res,dmin):" % (k, lin, ctl, nsol, N), [(float('%.2e' % b[0]), float('%.2e' % b[1])) for b in best[:5]])
