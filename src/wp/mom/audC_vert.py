# audC_vert.py -- vertical case (all cubic members share a = a0) of collinear + 2t+1, residual all-B subcase.
# Members (a0, c), c in the B-slice (max(4, a0 - sqrt(8 a0)), min(a0 - 4, sqrt(8 a0))), a0 in (8, 32).
# INT(t,v) holds iff c_t - c_b <= F_t(c_t) = (4a + c_t(a - c_t))/(2(a - c_t + 2));
# INT(b,w) holds iff c_t - c_b <= F_b(c_b) = (4a + c_b(a - c_b))/(2(c_b + 2)).
# Check: G(a, c_t) := F_t(c_t) - (c_t - c_lo(a)) >= 0 on the slice (then INT(t,v) alone suffices). Grid + exact-ish minimum.
import numpy as np, mpmath as mp
worst = (np.inf, None)
for a in np.linspace(8.0001, 31.9999, 4001):
    lo = max(4.0, a - np.sqrt(8*a)); hi = min(a - 4, np.sqrt(8*a))
    if hi <= lo: continue
    ct = np.linspace(lo, hi, 2001)
    G = (4*a + ct*(a - ct))/(2*(a - ct + 2)) - (ct - lo)
    k = np.argmin(G)
    if G[k] < worst[0]: worst = (G[k], (a, ct[k], lo, hi))
print("min over vertical B-slices of F_t(c_t) - (c_t - c_lo) =", worst)
# direct check with moduli for random vertical all-B pairs
import random
rnd = random.Random(3); bad = 0; n = 0
def f(A, C, t): return 2*t + 1 + t*(t+1)*(A*t + C)
def ur(A, B, C):
    d = mp.sqrt(B*B - 4*A*C); r = (-B + d)/(2*A); return r if mp.im(r) > 0 else (-B - d)/(2*A)
mp.mp.dps = 30
for _ in range(20000):
    a = mp.mpf(rnd.uniform(8, 32)); lo = max(4, a - mp.sqrt(8*a)); hi = min(a - 4, mp.sqrt(8*a))
    if hi <= lo: continue
    cb, ct = sorted([lo + (hi - lo)*mp.mpf(rnd.random()) for _ in range(2)])
    n += 1; v = ur(a, ct, 2)
    if not (abs(f(a, cb, v)) <= 1 and abs(2*v + 1) < 1): bad += 1
print("random vertical all-B pairs: %d, INT(t,v) failures: %d" % (n, bad))
