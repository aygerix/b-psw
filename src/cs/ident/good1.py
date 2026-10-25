#!/usr/bin/env python3
r"""src/cs/ident/good1.py -- count Korselt-type universal sets S having AT LEAST ONE good member
(N_d(S) = 0, i.e. member d is Lucas-automatic (F-form) / Fermat-automatic (L-form) for every m).
For a good d every e in S \ {d} lies in comp[d] (radical condition), so enumerate (k-1)-subsets of comp[d].
Usage: good1.py FORM Dmax kmax [maxcomp]"""
import sys
from math import gcd
from functools import reduce
from itertools import combinations
form, Dmax, kmax = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
maxcomp = int(sys.argv[4]) if len(sys.argv) > 4 else 60
fac = (lambda d, e: d - 2*e) if form == 'F' else (lambda d, e: 2*e - d)
def rad_divides(x, d):
    x = abs(x)
    if x == 0: return False
    g = gcd(x, d)
    while g > 1:
        while x % g == 0: x //= g
        g = gcd(x, d)
    return x == 1
lcm = lambda S: reduce(lambda a, b: a*b//gcd(a, b), S)
tot = {k: 0 for k in range(3, kmax+1)}; kor = {k: [] for k in range(3, kmax+1)}; skipped = 0; maxc = 0
for d in range(2, Dmax+1):
    C = [e for e in range(1, 2*Dmax+1) if e != d and rad_divides(fac(d, e), d)]
    maxc = max(maxc, len(C))
    if len(C) > maxcomp: skipped += 1; continue
    for k in range(3, kmax+1):
        if form == 'L' and k % 2 == 0: continue
        for R in combinations(C, k-1):
            p = 1
            for e in R: p *= fac(d, e)
            if p != d**(k-1): continue
            S = sorted(R + (d,))
            if reduce(gcd, S) != 1: continue
            tot[k] += 1
            if sum(S) % lcm(S) == 0: kor[k].append(S)
print(f"form {form} d<= {Dmax} (others <= {2*Dmax}) kmax {kmax}; max |comp| = {maxc}; d skipped (|comp|>{maxcomp}) = {skipped}")
for k in range(3, kmax+1):
    print(f"  k={k}: sets with >=1 good member (gcd 1) = {tot[k]}; of these Korselt-type (lcm | sum) = {len(kor[k])}: {kor[k][:10]}")
