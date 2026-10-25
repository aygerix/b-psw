#!/usr/bin/env python3
"""Control for oc_mersprod.py Part A: brute force over all subsets (|S| >= 2) of the first 12 odd Mersenne primes
(q <= 521): build n and run the actual strong base-2 test from src/bpsw.py; compare with the q_i | n - 1 criterion.
Also a positive control of the criterion on random sets of primes p with odd ord_p 2."""
import sys, os, itertools, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from bpsw import strong_fermat2
Q = [3,5,7,13,17,19,31,61,89,107,127,521]
agree = dis = spsp = 0
for k in range(2, len(Q) + 1):
    for S in itertools.combinations(Q, k):
        n = 1
        for q in S: n *= (1 << q) - 1
        crit = all((n - 1) % q == 0 for q in S)
        sf = strong_fermat2(n)
        spsp += sf
        if crit == sf: agree += 1
        else: dis += 1; print('DISAGREE', S, crit, sf)
print(f'Mersenne subsets q<=521: {agree + dis} checked, criterion agrees {agree}, disagrees {dis}, strong base-2 psp {spsp}')
# positive control: primes with ord_p 2 odd, products in class j1 = 0
from sympy import isprime, n_order
pool = [p for p in range(3, 20000) if isprime(p) and n_order(2, p) % 2 == 1]
hits = 0; agree2 = 0; tot = 0
random.seed(1)
for _ in range(200000):
    S = random.sample(pool, 3); n = S[0]*S[1]*S[2]
    crit = all((n - 1) % n_order(2, p) == 0 for p in S)
    sf = strong_fermat2(n); tot += 1; agree2 += (crit == sf); hits += sf
print(f'positive control: {tot} random triples from {len(pool)} primes with odd ord: criterion agrees {agree2}/{tot}, spsp hits {hits}')
# exhaustive pairs (has true positives, e.g. 2047 = 23*89)
od = {p: n_order(2, p) for p in pool}
tot = agree3 = hits = 0; ex = []
for a, b in itertools.combinations(pool, 2):
    n = a * b; crit = (n - 1) % od[a] == 0 and (n - 1) % od[b] == 0
    sf = strong_fermat2(n); tot += 1; agree3 += (crit == sf); hits += sf
    if sf and len(ex) < 5: ex.append((a, b))
print(f'pair control: {tot} pairs, criterion agrees {agree3}/{tot}, spsp hits {hits}, e.g. {ex}')
