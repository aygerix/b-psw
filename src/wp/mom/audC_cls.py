# audC_cls.py -- which robust endpoint test passes, tabulated by the regions of the two endpoints.
# i = endpoint with smaller a, j = endpoint with larger a. Regions (O = linear member vs the ellipses of X):
#   N: c <= 4 and a - c <= 4 (O outside both);  W: a - c > 4, c <= 4 (O inside E^w);  V: c > 4, a - c <= 4;  B: both.
# Tests as in audC_rob.py. Also records the states of the other endpoint relative to each (X, z):
#   H (E >= 0: closed outer half-plane), I (|factor| <= 1), X (outside ellipse but E < 0).
import sys, numpy as np, collections
rng = np.random.default_rng(int(sys.argv[1]))
def fv(a, c, t): return 2*t + 1 + t*(t+1)*(a*t + c)
def up(A, B, C):
    r = np.roots([A, B, C]); return r[np.argmax(r.imag)]
def inD(a, c): return a > 0 and c*c < 8*a and (a-c)**2 < 8*a
def region(a, c): return "NWVB"[(a - c > 4) + 2*(c > 4)]
def state(a, c, other, zn):
    w = up(a, a + c, c + 2); v = up(a, c, 2)
    F = (lambda Q: fv(Q[0], Q[1], w)) if zn == "w" else (lambda Q: -fv(Q[0], Q[1], v))
    val = F(other); E = val.real - 1
    st = "H" if E >= 0 else ("I" if abs(val) <= 1 else "X")
    Oin = abs(F((0.0, 0.0))) < 1
    ok = (st == "H" and not Oin) or (st == "I" and Oin)
    return st, Oin, ok
tab = collections.defaultdict(lambda: collections.Counter()); tot = collections.Counter(); unc = 0
T = int(sys.argv[2])
def samp():
    while True:
        a = 32*rng.random()**rng.choice([1, 2, 0.3]); lo, hi = a - np.sqrt(8*a), np.sqrt(8*a)
        c = lo + (hi - lo)*rng.random()
        if inD(a, c): return (a, c)
for _ in range(T):
    Q1, Q2 = samp(), samp()
    if rng.random() < 0.3:   # short segments too
        s = rng.choice([1e-3, 1e-2, 1e-1]); Q2 = (Q1[0] + s*(Q2[0]-Q1[0]), Q1[1] + s*(Q2[1]-Q1[1]))
    Qi, Qj = (Q1, Q2) if Q1[0] <= Q2[0] else (Q2, Q1)
    key = region(*Qi) + region(*Qj); tot[key] += 1
    passes = []
    for X, Y, tag in ((Qi, Qj, "i"), (Qj, Qi, "j")):
        for zn in "wv":
            st, Oin, ok = state(*X, Y, zn)
            if ok: passes.append(("EXT" if st == "H" else "INT") + tag + zn)
    if not passes: unc += 1
    for p in passes: tab[key][p] += 1
print("pairs %d, uncovered %d" % (T, unc))
for key in sorted(tot):
    n = tot[key]; always = [p for p, c in tab[key].items() if c == n]
    print("%s n=%6d  always-passing: %-30s  counts: %s" % (key, n, ",".join(sorted(always)), dict(sorted(tab[key].items()))))
