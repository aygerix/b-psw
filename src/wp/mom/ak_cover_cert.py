# ak_cover_cert.py -- verify a fully covered configuration in 60-digit arithmetic with RATIONAL member data.
# For each member j and each upper special point zeta (root of g_j or h_j), print min_{i != j} |f_i(zeta)| (< 1 = covered).
import sys, mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 60
A = [Fr(x) for x in sys.argv[1].split(',')]; C = [Fr(x) for x in sys.argv[2].split(',')]
n = len(A)
def f(a, c, t): return 2*t + 1 + t*(t+1)*(mp.mpf(a.numerator)/a.denominator*t + mp.mpf(c.numerator)/c.denominator)
worst = mp.mpf(0)
for j in range(n):
    a, c = A[j], C[j]
    dg = (a - c)**2 - 8*a; dh = c**2 - 8*a
    assert dg < 0 and dh < 0 and 0 < a < 32, "member %d not admissible" % j
    am, cm = mp.mpf(a.numerator)/a.denominator, mp.mpf(c.numerator)/c.denominator
    w = (-(am + cm) + mp.sqrt(mp.mpf(dg.numerator)/dg.denominator))/(2*am)
    v = (-cm + mp.sqrt(mp.mpf(dh.numerator)/dh.denominator))/(2*am)
    assert abs(f(a, c, w) - 1) < mp.mpf(10)**-50 and abs(f(a, c, v) + 1) < mp.mpf(10)**-50
    for name, z in (("w", w), ("v", v)):
        vals = [abs(f(A[i], C[i], z)) for i in range(n) if i != j]
        mn = min(vals); worst = max(worst, mn)
        print("member %d (a=%s, c=%s) point %s: min_{i!=j} |f_i| = %s" % (j, a, c, name, mp.nstr(mn, 12)))
print("max over special points of min coverage |f| = %s  (< 1: FULLY COVERED)" % mp.nstr(worst, 12))
