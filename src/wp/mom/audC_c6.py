# audC_c6.py -- evaluate the single-point straddle pattern at the line's C6 example (3 cubics, no linear member).
import mpmath as mp
mp.mp.dps = 50
a = [mp.mpf(2834289)/100000, mp.mpf(2836507)/100000, mp.mpf(1348957)/50000]
c = [mp.mpf(685267)/50000, mp.mpf(1449513)/100000, mp.mpf(335351)/25000]
def f(A, C, t): return 2*t + 1 + t*(t+1)*(A*t + C)
def ur(A, B, C):
    d = mp.sqrt(B*B - 4*A*C); r = (-B + d)/(2*A); return r if mp.im(r) > 0 else (-B - d)/(2*A)
for lin in (False, True):
    print("with linear member:", lin)
    for i in range(3):
        w = ur(a[i], a[i] + c[i], c[i] + 2); v = ur(a[i], c[i], 2)
        for nm, z in (("w", w), ("v", v)):
            mods = [abs(f(a[l], c[l], z)) for l in range(3) if l != i] + ([abs(2*z + 1)] if lin else [])
            st = "straddle" if (min(mods) < 1 < max(mods)) else ("all>=1" if min(mods) >= 1 else "all<=1")
            print("  member %d %s: moduli %s -> %s" % (i, nm, [mp.nstr(x, 8) for x in mods], st))
