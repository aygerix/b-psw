#!/usr/bin/env python3
"""src/cs/e3b/advchk/dpen_rat.py -- the real Sum m = 6 point on the reflection-symmetric DOUBLE-PENCIL Korselt component
(r_2 = r_3 = r, r_4 = r_5 = -r, f_4(x) = f_2(-x), f_5(x) = f_3(-x); excluded by the E3b nondeg filter).
Structure 1|1|2|34 (+ mirror): member 1 Lucas at +-w, members 3,4 share phi, members 2,5 share -phi.
Square system in (c, l2, l3, r, w, phi); polish to 80 digits; exact Lucas recount; rationality of c, l2 r, l3 r, r^2.
Usage: dpen_rat.py c l2 l3 r w phi   (starting values, real parts)"""
import sys, mpmath as mp
mp.mp.dps = 80
c0, l20, l30, r0, w0, p0 = [mp.mpc(eval(a)) for a in sys.argv[1:7]]

def fam(c, l2, l3, r):
    A = lambda x: x**2 - 1
    fs = [lambda x: 1 + A(x) * c,
          lambda x: 1 + A(x) * l2 * (x - r),
          lambda x: 1 + A(x) * l3 * (x - r),
          lambda x: 1 - A(x) * l2 * (x + r),
          lambda x: 1 - A(x) * l3 * (x + r)]
    return fs

def eqs(c, l2, l3, r, w, p):
    fs = fam(c, l2, l3, r)
    F = lambda x: fs[0](x) * fs[1](x) * fs[2](x) * fs[3](x) * fs[4](x)
    kor = (fs[0](r) * fs[3](r) * fs[4](r) - 1) / (r**2 - 1)
    return [kor, fs[0](w) + 1, F(w) + 1, fs[2](p) + 1, fs[3](p) + 1, F(p) + 1]

sol = mp.findroot(eqs, (c0, l20, l30, r0, w0, p0))
c, l2, l3, r, w, p = sol
print("residual", mp.nstr(mp.norm(mp.matrix(eqs(*sol))), 5))
print("imag parts", [mp.nstr(mp.im(v), 3) for v in sol])
# exact-ish Lucas recount at 80 digits: roots of f_j + 1 and F + 1 there
fs = fam(c, l2, l3, r)
F = lambda x: fs[0](x) * fs[1](x) * fs[2](x) * fs[3](x) * fs[4](x)
coefs = [[c, 0, 2 - c], [l2, -l2 * r, -l2, l2 * r + 2], [l3, -l3 * r, -l3, l3 * r + 2],
         [-l2, -l2 * r, l2, l2 * r + 2], [-l3, -l3 * r, l3, l3 * r + 2]]
gd = []
for j, cf in enumerate(coefs):
    rts = mp.polyroots(cf, maxsteps=200, extraprec=200)
    gd.append(sum(1 for x in rts if abs(F(x) + 1) < mp.mpf(10)**-60))
print("gcd degrees", gd, "sum", sum(gd))
# also check that f_j are pairwise distinct and the Korselt divisibilities hold numerically for all members
print("Korselt K at r (members 2,3) and -r (members 4,5):", mp.nstr(abs(fs[0](r) * fs[3](r) * fs[4](r) - 1), 3),
      mp.nstr(abs(fs[0](-r) * fs[1](-r) * fs[2](-r) - 1), 3))
for name, v in [('c', c), ('l2 r', l2 * r), ('l3 r', l3 * r), ('r^2', r**2)]:
    vr = mp.re(v)
    out = []
    for d in (1, 2, 3, 4):
        q = mp.findpoly(vr, d, maxcoeff=10**12)
        if q:
            out = (d, q); break
    print(f"{name} = {mp.nstr(vr, 30)}  minpoly(deg<=4, |coef|<=1e12): {out if out else 'none'}")
