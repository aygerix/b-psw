# audC_res.py -- residual cases of the collinear + 2t+1 extension after Lemma A (i in cl N) and Lemma B (i in W => EXTiv or j in W).
# i = smaller-a endpoint, j = larger-a endpoint. Regions as in audC_cls.py.
# Residual case WW: i, j in W and E^v_i(Q_j) < 0.   Residual case Bx: i in B (then a_j >= a_i > 8, j in W, V or B).
# For each residual pair record which of the 8 robust tests pass; report the tests that pass on ALL sampled pairs,
# and the minimum over pairs of the best margin among a candidate set (adversarial, Nelder-Mead).
import sys, numpy as np, collections
from scipy.optimize import minimize
rng = np.random.default_rng(int(sys.argv[1]))
def fv(a, c, t): return 2*t + 1 + t*(t+1)*(a*t + c)
def up(A, B, C):
    r = np.roots([A, B, C]); return r[np.argmax(r.imag)]
def inD(a, c): return 0 < a < 32 and c*c < 8*a and (a-c)**2 < 8*a
def region(a, c): return "NWVB"[(a - c > 4) + 2*(c > 4)]
def margins(Qi, Qj):
    out = {}
    for tag, (A, B) in (("i", (Qi, Qj)), ("j", (Qj, Qi))):
        a, c = A; w = up(a, a + c, c + 2); v = up(a, c, 2)
        for zn, F in (("w", lambda Q: fv(Q[0], Q[1], w)), ("v", lambda Q: -fv(Q[0], Q[1], v))):
            Z = F(B) - 1; FO = F((0.0, 0.0))
            out["EXT"+tag+zn] = min(Z.real/abs(Z), np.log(abs(FO)))
            out["INT"+tag+zn] = min(-np.log(abs(1 + Z)), -np.log(abs(FO)))
    return out
def Ev(Qi, Qj):
    a, c = Qi; v = up(a, c, 2); return (-fv(Qj[0], Qj[1], v)).real - 1
def samp(reg):
    while True:
        a = 8 + 24*rng.random() if reg == "B" else 32*rng.random()
        lo, hi = a - np.sqrt(8*a), np.sqrt(8*a); c = lo + (hi - lo)*rng.random()
        if inD(a, c) and (reg is None or region(a, c) == reg): return (a, c)
T = int(sys.argv[2]); cases = {"WW": collections.Counter(), "Bx": collections.Counter()}; n = collections.Counter(); unc = collections.Counter()
for _ in range(T):
    # WW residual
    Q1, Q2 = samp("W"), samp("W")
    if rng.random() < 0.3: s = rng.choice([1e-3, 1e-2, 1e-1]); Q2 = (Q1[0] + s*(Q2[0]-Q1[0]), Q1[1] + s*(Q2[1]-Q1[1]))
    Qi, Qj = (Q1, Q2) if Q1[0] <= Q2[0] else (Q2, Q1)
    if region(*Qi) == "W" and region(*Qj) == "W" and Ev(Qi, Qj) < 0:
        n["WW"] += 1; m = margins(Qi, Qj)
        for k_, v_ in m.items():
            if v_ >= 0: cases["WW"][k_] += 1
        if max(m.values()) < 0: unc["WW"] += 1
    # B residual: i in B, j anywhere in D with a_j >= a_i
    Qi = samp("B"); Qj = samp(None)
    if rng.random() < 0.3: s = rng.choice([1e-3, 1e-2, 1e-1]); Qj = (Qi[0] + s*(Qj[0]-Qi[0]), Qi[1] + s*(Qj[1]-Qi[1]))
    if Qj[0] >= Qi[0] and inD(*Qj):
        key = "B" + region(*Qj); n[key] += 1; cases.setdefault(key, collections.Counter()); m = margins(Qi, Qj)
        for k_, v_ in m.items():
            if v_ >= 0: cases[key][k_] += 1
        if max(m.values()) < 0: unc[key] += 1
for key in sorted(n):
    always = [t for t, c in cases[key].items() if c == n[key]]
    print("%s: n=%d uncovered=%d always=%s counts=%s" % (key, n[key], unc[key], always, dict(sorted(cases[key].items()))))
