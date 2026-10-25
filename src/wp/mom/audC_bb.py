# audC_bb.py -- structure of the BB residual (i, j in B, a_i <= a_j): split by the slope s = (c_j - c_i)/(a_j - a_i).
#   s <= 0 (c_j <= c_i): test INT(j,w);  s >= 1 (a_j - c_j <= a_i - c_i): test INT(j,v);  0 < s < 1: which tests pass.
# Same for BW (i in B, j in W) and WW-residual, with the slope.
import sys, numpy as np, collections
rng = np.random.default_rng(int(sys.argv[1]))
def fv(a, c, t): return 2*t + 1 + t*(t+1)*(a*t + c)
def up(A, B, C):
    r = np.roots([A, B, C]); return r[np.argmax(r.imag)]
def inD(a, c): return 0 < a < 32 and c*c < 8*a and (a-c)**2 < 8*a
def region(a, c): return "NWVB"[(a - c > 4) + 2*(c > 4)]
def sampB():
    while True:
        a = 8 + 24*rng.random(); lo, hi = max(4, a - np.sqrt(8*a)), min(a - 4, np.sqrt(8*a))
        if hi > lo: return (a, lo + (hi - lo)*rng.random())
def tests(Qi, Qj):
    out = {}
    for tag, (A, B) in (("i", (Qi, Qj)), ("j", (Qj, Qi))):
        a, c = A; w = up(a, a + c, c + 2); v = up(a, c, 2)
        for zn, F in (("w", lambda Q: fv(Q[0], Q[1], w)), ("v", lambda Q: -fv(Q[0], Q[1], v))):
            Z = F(B) - 1; FO = F((0.0, 0.0))
            if Z.real >= 0 and abs(FO) >= 1: out["EXT"+tag+zn] = 1
            if abs(1 + Z) <= 1 and abs(FO) <= 1: out["INT"+tag+zn] = 1
    return out
cnt = collections.defaultdict(collections.Counter); n = collections.Counter()
for _ in range(int(sys.argv[2])):
    Q1, Q2 = sampB(), sampB()
    if rng.random() < 0.3: s = rng.choice([1e-3, 1e-2, 1e-1]); Q2 = (Q1[0] + s*(Q2[0]-Q1[0]), Q1[1] + s*(Q2[1]-Q1[1]))
    Qi, Qj = (Q1, Q2) if Q1[0] <= Q2[0] else (Q2, Q1)
    if Qj[0] == Qi[0]: continue
    s = (Qj[1] - Qi[1])/(Qj[0] - Qi[0]); key = "s<=0" if s <= 0 else ("s>=1" if s >= 1 else "0<s<1")
    n[key] += 1
    for t in tests(Qi, Qj): cnt[key][t] += 1
for key in n:
    print("BB %s: n=%d always=%s counts=%s" % (key, n[key], [t for t, c in cnt[key].items() if c == n[key]], dict(sorted(cnt[key].items()))))
