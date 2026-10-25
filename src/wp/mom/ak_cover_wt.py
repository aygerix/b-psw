# ak_cover_wt.py -- for a configuration, per member P list the vectors (log|f_Q(w_P)|, log|f_Q(v_P)|), Q != P, and decide
# whether some NONNEGATIVE weights (lam, mu) make all lam*x + mu*y > 0 (a weighted single-member extremal argument at P),
# and whether ANY real weights do (i.e. 0 not in the convex hull).   60-digit arithmetic, rational input.
import sys, mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 60
A = [Fr(x) for x in sys.argv[1].split(',')]; C = [Fr(x) for x in sys.argv[2].split(',')]; n = len(A)
q = lambda r: mp.mpf(r.numerator)/r.denominator
def f(i, t): return 2*t + 1 + t*(t+1)*(q(A[i])*t + q(C[i]))
for j in range(n):
    a, c = q(A[j]), q(C[j])
    w = (-(a + c) + mp.sqrt((a - c)**2 - 8*a))/(2*a); v = (-c + mp.sqrt(c**2 - 8*a))/(2*a)
    vecs = [(mp.log(abs(f(i, w))), mp.log(abs(f(i, v)))) for i in range(n) if i != j]
    # nonneg weights: need a direction theta in [0, pi/2] with cos*x + sin*y > 0 for all vectors; scan finely
    best = max(min(mp.cos(th)*x + mp.sin(th)*y for (x, y) in vecs) for th in [mp.pi/2*s/2000 for s in range(2001)])
    bestany = max(min(mp.cos(th)*x + mp.sin(th)*y for (x, y) in vecs) for th in [2*mp.pi*s/8000 for s in range(8000)])
    print("member %d: vectors %s ; best nonneg-weight margin %s ; best any-weight margin %s" % (j, [(mp.nstr(x, 6), mp.nstr(y, 6)) for x, y in vecs], mp.nstr(best, 6), mp.nstr(bestany, 6)))
