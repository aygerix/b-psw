#!/usr/bin/env python3
"""src/cs/ident/two_bad.py -- search for universal forms with EXACTLY two 'bad' members.
For a set S (gcd 1) of cofactors, member d is GOOD iff its fixed residual N_d(S) vanishes:
  F (Chernick, p = 1 + d m):        prod_{e!=d} (d - 2e)        = d^(k-1)
  L (Lucas-Carmichael, p = d m - 1): prod_{e!=d} (2e - d)        = d^(k-1)
  T (p = 2 (d m)^2 - 1):             prod_{e!=d} (2e^2 - d^2)   = d^(2(k-1))
Sign lemma (PROVED): min(S) and max(S) are always bad.  A good member d forces every factor to be nonzero with
radical dividing rad(d) (the product is +- a power of d).  We search cliques of mutually compatible intermediate
members, then the extremes, then test exact equality and the additive (Korselt-type) condition.
Usage: two_bad.py FORM Dmax kmax"""
import sys
from math import gcd
from functools import reduce
form, Dmax, kmax = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])

def fac(d, e):
    if form == 'F': return d - 2*e
    if form == 'L': return 2*e - d
    return 2*e*e - d*d

def rad_divides(x, d):
    x = abs(x)
    if x == 0: return False
    g = gcd(x, d)
    while g > 1:
        while x % g == 0: x //= g
        g = gcd(x, d)
    return x == 1

comp = [set() for _ in range(Dmax+1)]   # comp[d] = {e : factor(d,e) nonzero with rad | rad(d)}
for d in range(1, Dmax+1):
    for e in range(1, Dmax+1):
        if e != d and rad_divides(fac(d, e), d):
            comp[d].add(e)
power = (lambda k: k-1) if form != 'T' else (lambda k: 2*(k-1))

def good(S, d):
    k = len(S); p = 1
    for e in S:
        if e != d: p *= fac(d, e)
    return p == d**power(k)

def additive_ok(S):
    L = reduce(lambda a, b: a*b//gcd(a, b), S)
    if form == 'T': return sum(x*x for x in S) % (L*L) == 0
    return sum(S) % L == 0

found = {k: [] for k in range(3, kmax+1)}
cnt_exact = {k: 0 for k in range(3, kmax+1)}
def extend(I, cand):
    # I: sorted clique of intermediates (mutually compatible); cand: candidates > max(I) compatible with all of I both ways
    if len(I) >= 1:
        lo = min(I); hi = max(I)
        A = set.intersection(*[comp[d] for d in I])
        for a in [x for x in A if x < lo]:
            for b in [x for x in A if x > hi]:
                S = sorted([a] + I + [b]); k = len(S)
                if k > kmax: continue
                if reduce(gcd, S) != 1: continue
                if form != 'F' and k % 2 == 0: continue
                if all(good(S, d) for d in I):
                    cnt_exact[k] += 1
                    if additive_ok(S): found[k].append(S)
    if len(I) + 2 >= kmax: return
    for e in sorted(cand):
        newc = {x for x in cand if x > e and x in comp[e] and e in comp[x]}
        extend(I + [e], newc)

for d in range(1, Dmax+1):
    extend([d], {x for x in comp[d] if x > d and d in comp[x]})
print(f"form {form} Dmax {Dmax} kmax {kmax}")
for k in range(3, kmax+1):
    print(f"  k={k}: exactly-two-bad sets (gcd 1) = {cnt_exact[k]}; with additive condition = {len(found[k])}: {found[k][:12]}")
