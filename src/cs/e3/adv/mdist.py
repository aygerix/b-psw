# Moduli M = lcm of the content denominators per (rho, d) in the |a|, b <= 30, d in {1,2} grid; identifies the
# (rho, d) pairs that rho_scan.gp skipped (M > MMAX) -- the 8 families with rows in v5_scan but none in rho_scan.
import sys, math, os
HERE = os.path.dirname(os.path.abspath(__file__))
from fractions import Fraction as Fr
sys.argv = ['x', os.path.join(HERE, '../../..'), '1', '1', '1', '1', '/dev/null', '--estimate']
exec(open(os.path.join(HERE, 'v5_scan.py')).read().split('t0 = time.time()')[0])
skipped8 = {(-13,25),(-16,25),(-17,27),(-22,21),(-23,11),(-26,21),(26,17),(29,15)}
Ms = []
for b in range(1, 31):
    for a in range(-30, 31):
        if a == 0 or math.gcd(a, b) != 1 or 2 * a == b: continue
        f = family(Fr(a, b))
        for d in (1, 2):
            P = [intpoly(fj, d) for fj in f]
            M = 1
            for _, L in P: M = M * L // math.gcd(M, L)
            Ms.append((M, a, b, d, (a, b) in skipped8 and d == 2))
Ms.sort()
print("families", len(Ms), "max M", Ms[-1])
print("min M among the 8 missing (d=2):", min(m for m in Ms if m[4]))
for thr in (10**5, 2*10**5, 5*10**5, 10**6, 2*10**6, 5*10**6):
    print(thr, "families with M >", thr, ":", sum(1 for m in Ms if m[0] > thr), " of which in the missing 8:", sum(1 for m in Ms if m[0] > thr and m[4]))
print([m for m in Ms if m[4]])
