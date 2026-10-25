#!/usr/bin/env python3
"""src/cs/e3b/fsp_mp.py -- E3b: high-precision refinement (mpmath Newton, square FSP system) of FSP + 1 points and
the field test: the family is defined over Q iff T = t^2 and u_i = (t - r_i)/t are all rational (A = x^2 - d, scaling
x -> sqrt(d) x).  Prints the minimal-polynomial degree (mpmath.findpoly, deg <= 12, |coeff| <= 1e15) of T and of
each u_i, for every REAL point in the input (fsp_num.py output).  Usage: fsp_mp.py extra file [dps]"""
import sys, mpmath as mp, numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0] if '/' in __file__ else '.')
from fsp_num import FSP

def F(model, z):
    zz = np.array([z], dtype=object)
    return [mp.mpc(v) for v in model.F(zz)[0]]

def newton(model, z, dps):
    mp.mp.dps = dps
    z = [mp.mpc(v) for v in z]; n = len(z)
    for it in range(80):
        f = F(model, z); nf = mp.norm(mp.matrix(f))
        if nf < mp.mpf(10)**(-dps + 10):
            break
        h = mp.mpf(10)**(-dps // 2)
        J = mp.matrix(len(f), n)
        for i in range(n):
            zp = list(z); zp[i] += h; zm = list(z); zm[i] -= h
            fp = F(model, zp); fm = F(model, zm)
            for r in range(len(f)):
                J[r, i] = (fp[r] - fm[r]) / (2 * h)
        d = mp.lu_solve(J, -mp.matrix(f)) if len(f) == n else mp.lu_solve(J.H * J, -(J.H * mp.matrix(f)))
        z = [z[i] + d[i] for i in range(n)]
    return z, nf

def mindeg(v, maxdeg=12):
    v = mp.re(v)
    for d in range(1, maxdeg + 1):
        p = mp.findpoly(v, d, maxcoeff=10**15)
        if p:
            return d
    return '>%d' % maxdeg

extra = [e for e in sys.argv[1].split(',') if e]; fn = sys.argv[2]; dps = int(sys.argv[3]) if len(sys.argv) > 3 else 80
model = FSP(extra)
for line in open(fn):
    d, gd, zs = eval(line)
    if max(abs(complex(v).imag) for v in zs[:5]) > 1e-8:
        continue
    z, nf = newton(model, zs, dps)
    t = z[0]; rr = z[1:5]
    inv = [t**2] + [(t - r) / t for r in rr]
    degs = [mindeg(v) for v in inv]
    print(f"gcd {gd} resid {mp.nstr(nf, 3)}  T = {mp.nstr(mp.re(inv[0]), 25)}  min-poly degrees [T, u2..u5] = {degs}", flush=True)
