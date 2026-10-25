"""Reality-restricted sweep: for rational shapes on a grid (e1 > 0; e1 = 0 is the proved symmetric case), compute the reality
interval I (numerically, then an outward rational superset [LO, HI]), a constant alpha by LP over 7-member configurations in
[LO, HI]^7, and write a gp batch file for the EXACT branch-and-bound (pc7_cert_bb.gp).  Usage: python3 pc7_cert_rsweep.py step out.gp"""
import sys, numpy as np
from fractions import Fraction as Fr
from scipy.optimize import linprog
from pc7_cert_dlib import load, Evec
at = load(); rng = np.random.default_rng(23); step = Fr(sys.argv[1])
def interval(r1, r2):
    e1, e2 = r1 + r2, r1 * r2; Dl = (1 - r1**2) * (1 - r2**2); lam, mu = 2 * e1 / Dl, 2 * (1 - e1**2 + e2) / Dl
    q1 = np.roots([lam**2, -2 * lam * e1 - 4 * mu, e1**2 - 4 * e2]); q2 = np.roots([lam**2, -4 * (mu + lam * e1), 4])
    if np.any(np.abs(np.imag(q1)) > 1e-12) or np.any(np.abs(np.imag(q2)) > 1e-12): return None
    lo = max(np.min(np.real(q1)), np.min(np.real(q2))); hi = min(np.max(np.real(q1)), np.max(np.real(q2)))
    return (lo, hi) if lo < hi and hi > 0 else None
def outward(lo, hi):
    # rational LO <= lo and HI >= hi, within relative 1e-4
    L = Fr(lo * (1 - 1e-4)).limit_denominator(10**5); H = Fr(hi * (1 + 1e-4)).limit_denominator(10**5)
    while float(L) > lo: L -= Fr(1, 10**6)
    while float(H) < hi: H += Fr(1, 10**6)
    return L, H
out = open(sys.argv[2], 'w'); out.write('read("pc7_cert_bb.gp");\nRES = List();\n')
vals = [Fr(-95, 100) + k * step for k in range(int(Fr(190, 100) / step) + 1)]
for a in vals:
    for b in vals:
        if b < a or a + b <= 0: continue
        r1, r2 = float(a), float(b); I = interval(r1, r2)
        if I is None: print(f"r = ({r1:+.2f},{r2:+.2f}): reality interval EMPTY", flush=True); continue
        LO, HI = outward(*I)
        if LO <= 0: LO = Fr(1, 10**6)
        lo, hi = float(LO), float(HI); K, A, C = at(r1, r2)
        lg = lambda n: np.exp(rng.uniform(np.log(lo), np.log(hi), n))
        CF = [lg(7) for _ in range(1200)]
        for _ in range(1200):
            m = rng.integers(1, 4); v = lg(m); CF.append(v[rng.integers(0, m, 7)])
        g = np.exp(np.linspace(np.log(lo), np.log(hi), 9))
        for x in g:
            for y in g: CF.append(np.array([x] * 3 + [y] * 4)); CF.append(np.array([x] * 6 + [y])); CF.append(np.array([x] * 7))
        D = np.array([Evec(K, A, C, c) for c in CF]); ok = np.all(np.isfinite(D), axis=1); D = D[ok]; Dn = D / np.abs(D).sum(axis=1, keepdims=True)
        r = linprog(np.concatenate([np.zeros(8), [-1]]), A_ub=np.hstack([-Dn, np.ones((len(Dn), 1))]), b_ub=np.zeros(len(Dn)), bounds=[(-1, 1)] * 8 + [(None, None)], method='highs')
        al = r.x[:8]; t = r.x[-1]
        print(f"r = ({r1:+.2f},{r2:+.2f}): I = [{I[0]:.4g}, {I[1]:.4g}]  box [{LO}, {HI}]  LP margin {t:+.3e}", flush=True)
        if t <= 1e-9: continue
        alq = [Fr(x).limit_denominator(10**6) for x in al]
        out.write(f'SHAPE = [{a}, {b}]; LO = {LO}; HI = {HI}; ALPHA = [{", ".join(str(x) for x in alq)}];\n')
        out.write('OK = 0; foreach([6, 12, 24], n, N = n; if(bb() > 0, OK = n; break)); listput(RES, [SHAPE, OK]); print("RESULT ", SHAPE, " certified with N = ", OK);\n')
out.write('print("SUMMARY certified ", #select(x -> x[2] > 0, Vec(RES)), " of ", #RES);\n'); out.close()
