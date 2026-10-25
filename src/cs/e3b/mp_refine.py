#!/usr/bin/env python3
"""src/cs/e3b/mp_refine.py -- E3b: refine node-model solutions to high precision (mpmath Gauss-Newton) and test
whether the family is defined over Q.

Rationality test.  A shared-A family over Q has A = (x - h)^2 - d; the normalisation x -> sqrt(d) x + h gives
A = x^2 - 1 with  c' = d c,  l'_j = d^(3/2) l_j,  r'_j = (r_j - h)/sqrt(d).  Hence the normalised point comes from Q iff
  c', l'_j r'_j (= d l_j (r_j - h)) and r'_i r'_j (= (r_i-h)(r_j-h)/d) are all rational  (r'_j^2 in Q, ratios in Q).
We report the degree of the minimal polynomial (mpmath.findpoly, maxdeg 8, coefficients <= 1e12) of each invariant.
Usage: mp_refine.py STRUCT infile(lines from node_model --dump) [dps]
"""
import sys, mpmath as mp
sys.path.insert(0, __file__.rsplit('/', 1)[0] if '/' in __file__ else '.')
import numpy as np
from node_model import Node, parse
from refine import lucas_table

def Fmp(model, z):
    zz = np.array([z], dtype=object)
    with np.errstate(all='ignore'):
        out = model.F(zz)[0]
    return [mp.mpc(v) for v in out]

def newton(model, z, dps=60, iters=60):
    mp.mp.dps = dps
    z = [mp.mpc(v) for v in z]
    n = len(z)
    for it in range(iters):
        f = Fmp(model, z)
        nf = mp.norm(mp.matrix(f))
        if nf < mp.mpf(10)**(-dps + 8):
            break
        h = mp.mpf(10)**(-dps // 2)
        Jc = []
        for i in range(n):
            zp = list(z); zp[i] += h
            zm = list(z); zm[i] -= h
            fp = Fmp(model, zp); fm = Fmp(model, zm)
            Jc.append([(a - b) / (2 * h) for a, b in zip(fp, fm)])
        J = mp.matrix(len(f), n)
        for i in range(n):
            for r in range(len(f)):
                J[r, i] = Jc[i][r]
        JH = J.H
        d = mp.lu_solve(JH * J, -(JH * mp.matrix(f)))
        z = [z[i] + d[i] for i in range(n)]
    return z, nf

def invariants(model, z):
    zz = np.array([z], dtype=object)
    c = model.gval(zz, 1, mp.mpf(0))[0]
    lam = []; rr = [z[j] for j in range(4)]
    for j in range(2, 6):
        g0 = model.gval(zz, j, mp.mpf(0))[0]; g1 = model.gval(zz, j, mp.mpf(1))[0]
        lam.append(g1 - g0)
    inv = {'c': c}
    for j in range(4):
        inv['l%dr%d' % (j + 2, j + 2)] = lam[j] * rr[j]
    for i in range(4):
        for j in range(i, 4):
            inv['r%dr%d' % (i + 2, j + 2)] = rr[i] * rr[j]
    return inv, c, lam, rr

def mindeg(v, maxdeg=8):
    if abs(mp.im(v)) > mp.mpf(10)**(-30):
        return 'cplx'
    v = mp.re(v)
    for d in range(1, maxdeg + 1):
        p = mp.findpoly(v, d, maxcoeff=10**12)
        if p:
            return d if d > 1 else ('Q', mp.identify(v))
    return '>%d' % maxdeg

def main():
    struct = sys.argv[1]; fn = sys.argv[2]; dps = int(sys.argv[3]) if len(sys.argv) > 3 else 60
    model = Node(parse(struct))
    for line in open(fn):
        if not line.strip():
            continue
        st, gd, clr, zs = eval(line)
        if st != struct:
            continue
        z, nf = newton(model, zs, dps=dps)
        inv, c, lam, rr = invariants(model, z)
        allreal = all(abs(mp.im(v)) < mp.mpf(10)**(-30) for v in z)
        degs = {k: mindeg(v) for k, v in inv.items()}
        rational = all(isinstance(d, tuple) for d in degs.values())
        print(f"{struct} gcd {gd} residual {mp.nstr(nf, 3)} real={allreal} rational_family={rational}")
        print("   invariant min-poly degrees:", degs)
        print("   c =", mp.nstr(c, 20), " r =", [mp.nstr(r, 15) for r in rr])

if __name__ == "__main__":
    main()
