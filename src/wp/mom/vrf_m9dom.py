"""vrf line: independent exact check (Python Fractions, rational enclosure of sqrt D) of the domain inequalities
for the claimed k = 9 point (pair choice [1,3] of m9_symexact).  Domain for member (u, v):
  v > 1/32,  u^2 < 8v,  (1-u)^2 < 8v   [g, h without real roots]  <=>  (|u - 1/2| + 1/2)^2 < 8v.
With u = 1/2 +- sqrt(y): need y >= 0 and sqrt(y) + 1/2 < sqrt(8v).  We bound sqrt(y) by rational intervals
(no floating point) and check (sqrt(y)_hi + 1/2)^2 < 8 v, plus distinctness of the 9 members."""
from fractions import Fraction as Fr
from math import isqrt

D = Fr(40873015522726052437652941104413136372114406461, 426655372299127257812500000000000000000000000)
N = Fr(-3908730585275433747927744939801, 12548687420562566406250000000)
Dn = Fr(-286107999639694964577, 2362105867400012500)
c0 = Fr(3640457795177541991529, 138017997877500000000)
c1 = Fr(-3585828512, 360831367)
PV = [Fr(15129, 25000), Fr(671, 12500), Fr(6201, 100000), Fr(1193, 20000)]
vc = Fr(80537, 50000)
w2, w4 = Fr(2863, 10**5), Fr(1143, 10**4)

SC = 10**80


def sqrt_enclose(x):
    """rational lo <= sqrt(x) <= hi, exactly certified"""
    assert x >= 0
    n = isqrt((x.numerator * SC * SC) // x.denominator)
    lo, hi = Fr(n, SC), Fr(n + 1, SC)
    assert lo * lo <= x <= hi * hi
    return lo, hi


# D is not a rational square: numerator, denominator coprime -> both must be squares
num, den = D.numerator, D.denominator
print("D rational square?", isqrt(num) ** 2 == num and isqrt(den) ** 2 == den)
sDlo, sDhi = sqrt_enclose(D)
# y1 = (N - sqrt D)/Dn, Dn < 0 -> increasing in sqrt D
assert Dn < 0
y1lo, y1hi = (N - sDlo) / Dn, (N - sDhi) / Dn
y1lo, y1hi = min(y1lo, y1hi), max(y1lo, y1hi)
# y3 = c0 + c1 y1, c1 < 0
y3lo, y3hi = sorted([c0 + c1 * y1lo, c0 + c1 * y1hi])
print("y1 in [%.25f, %.25f], width %.3e" % (float(y1lo), float(y1hi), float(y1hi - y1lo)))
print("y3 in [%.25f, %.25f], width %.3e" % (float(y3lo), float(y3hi), float(y3hi - y3lo)))
# minimal polynomial residual sanity: (Dn y - N)^2 - D at the enclosure endpoints changes sign?
f = lambda y: (Dn * y - N) ** 2 - D
print("min poly sign change across y1 enclosure:", (f(y1lo) > 0) != (f(y1hi) > 0))

ok = True
for name, (ylo, yhi), v in [("y1", (y1lo, y1hi), PV[0]), ("y2", (w2 * w2, w2 * w2), PV[1]),
                            ("y3", (y3lo, y3hi), PV[2]), ("y4", (w4 * w4, w4 * w4), PV[3])]:
    if ylo <= 0:
        ok = False
        print(name, "not certified positive")
        continue
    slo, _ = sqrt_enclose(ylo)
    _, shi = sqrt_enclose(yhi)
    margin = 8 * v - (shi + Fr(1, 2)) ** 2
    print(name, " y > 0 certified; v > 1/32:", v > Fr(1, 32), "; certified lower bound of 8v - (sqrt y + 1/2)^2 = %.12f" % float(margin))
    ok &= margin > 0 and v > Fr(1, 32)
print("centre: 8 vc - 1/4 =", float(8 * vc - Fr(1, 4)), vc > Fr(1, 32))
ok &= 8 * vc > Fr(1, 4)
# distinctness: pair v's distinct and != vc; each pair has w != 0
print("pair v distinct and != vc:", len(set(PV + [vc])) == 5)
ok &= len(set(PV + [vc])) == 5
# c_i^2 < 8 a_i and (a_i - c_i)^2 < 8 a_i directly in (a, c) for the rational members (sanity of the reformulation)
for v, w in [(PV[1], w2), (PV[3], w4)]:
    for s in (1, -1):
        u = Fr(1, 2) + s * w
        a, c = 1 / v, u / v
        ok &= c * c < 8 * a and (a - c) ** 2 < 8 * a
        assert c * c < 8 * a and (a - c) ** 2 < 8 * a
print("ALL DOMAIN INEQUALITIES CERTIFIED:", ok)
