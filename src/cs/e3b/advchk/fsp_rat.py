#!/usr/bin/env python3
import os
"""src/cs/e3b/advchk/fsp_rat.py -- field test for the REAL FSP points of fsp_tu.py output files (independent code).
Polishes each real point to dps digits with mpmath Newton on the square fsp_tu system (object arrays), then reports for
T and each u_i the smallest-degree integer polynomial found by mpmath.findpoly (deg <= maxdeg, |coeff| <= maxcoeff).
A family is defined over Q iff T and all u_i are rational (Lemma B).  Usage: fsp_rat.py TYPE file [dps]"""
import sys, numpy as np, mpmath as mp
sys.path.insert(0, os.path.expanduser('~/b-psw/src/cs/e3b/advchk'))
import fsp_tu

typ = sys.argv[1]; fn = sys.argv[2]; dps = int(sys.argv[3]) if len(sys.argv) > 3 else 60
mp.mp.dps = dps

def Fv(z):
    Z = np.array([z], dtype=object)
    return [mp.mpc(v) for v in fsp_tu.system(Z, typ)[0]]

def newton(z):
    z = [mp.mpc(v) for v in z]; n = len(z)
    for it in range(60):
        f = Fv(z); nf = mp.norm(mp.matrix(f))
        if nf < mp.mpf(10) ** (-dps + 8):
            break
        h = mp.mpf(10) ** (-dps // 2)
        J = mp.matrix(n, n)
        for i in range(n):
            zp = list(z); zp[i] += h; zm = list(z); zm[i] -= h
            fp = Fv(zp); fm = Fv(zm)
            for r in range(n):
                J[r, i] = (fp[r] - fm[r]) / (2 * h)
        d = mp.lu_solve(J, -mp.matrix(f))
        z = [z[i] + d[i] for i in range(n)]
    return z, nf

def mindeg(v, maxdeg=12, maxcoeff=10**15):
    for d in range(1, maxdeg + 1):
        p = mp.findpoly(v, d, maxcoeff=maxcoeff)
        if p:
            return d, p
    return None

seen = 0
for line in open(fn):
    gd, z, near, sv = eval(line)
    z = np.array(z)
    if np.max(np.abs(z[:5].imag)) > 1e-8:
        continue
    seen += 1
    zz, nf = newton(list(z))
    vals = [mp.re(zz[0])] + [mp.re(w) for w in zz[1:5]]
    imag = max(abs(mp.im(w)) for w in zz[:5])
    degs = []
    for v in vals:
        r = mindeg(v)
        degs.append(r[0] if r else '>12')
    print(f"real pt {seen}: gd {gd} resid {mp.nstr(nf, 3)} imag {mp.nstr(imag, 3)} T = {mp.nstr(vals[0], 15)}  minpoly degrees (T, u2..u5): {degs}", flush=True)
print("real points tested:", seen)
