#!/usr/bin/env python3
"""Products of distinct Mersenne primes as BPSW candidates (Line OC, family 1). Exact.
Part A: every subset (|S| >= 2) of the 19 odd Mersenne exponents 3 <= q <= 4423:
  Fermat: all j1 = 0 (ord = q odd); need q_i | n - 1 for all i.
  Selfridge D from per-member Jacobi symbols; Lucas class j2 equal for all members (table oc_j2table.out),
  v2 room; survivors -> src/bpsw.py on the actual n.
Part B: Fermat-side census over all 51 known Mersenne exponents, |S| = 2..5 (n mod q_i only)."""
import sys, os, itertools, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from bpsw import bpsw_pass, jacobi
ALLQ = [2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,
        44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,
        24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933,136279841]
Q19 = [3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423]
J = {}
for line in open(sys.argv[1]):
    D, q, e, j = map(int, line.split()); J[(D, q)] = (e, j)
Dseq = []
D = 5
while abs(D) <= 103:
    Dseq.append(D); D = -D - 2 if D > 0 else -D + 2
def jacobi_char(D, r):
    # (D / M) for M == r mod 4|D|, r odd positive (Jacobi symbol is a character mod 4|D| in the bottom)
    return jacobi(D % r, r) if r > 1 else 1
t0 = time.time()
# Part A
mod_q = {qi: {qj: (pow(2, qj, qi) - 1) % qi for qj in Q19} for qi in Q19}
fermat_ok = 0; sel_ok = 0; class_ok = 0; bp = 0; nodd = 0
for k in range(2, len(Q19) + 1):
    for S in itertools.combinations(Q19, k):
        good = True
        for qi in S:
            r = 1
            for qj in S: r = r * mod_q[qi][qj] % qi
            if r != 1 % qi: good = False; break
        if not good: continue
        fermat_ok += 1
        # Selfridge D
        Dsel = None
        for D in Dseq:
            if int(abs(D) ** 0.5 + 0.5) ** 2 == abs(D) and D > 0: continue   # (square/n) never -1
            s = 1
            for q in S:
                e = J[(D, q)][0]
                s *= e
            if s == -1: Dsel = D; break
        if Dsel is None: print('  no D <= 103 for', S); continue
        sel_ok += 1
        js = set(J[(Dsel, q)][1] for q in S)
        if len(js) != 1: continue
        class_ok += 1
        n = 1
        for q in S: n *= (1 << q) - 1
        ok = bpsw_pass(n); bp += ok
        print('  CLASS PASS', S, 'D', Dsel, 'j2', js, 'bpsw', ok)
print(f'Part A: subsets {2**len(Q19) - 1 - len(Q19)}; Fermat (q_i | n-1) {fermat_ok}; Selfridge D found {sel_ok}; one Lucas class {class_ok}; BPSW {bp}; {time.time()-t0:.1f}s')
# Part B: Fermat-side census, all 51 odd exponents, |S| <= 5
QB = [q for q in ALLQ if q > 2]
modB = {qi: {qj: (pow(2, qj, qi) - 1) % qi for qj in QB} for qi in QB}
cnt = {}
surv = []
for k in range(2, 6):
    c = 0
    for S in itertools.combinations(QB, k):
        good = True
        for qi in S:
            r = 1
            for qj in S: r = r * modB[qi][qj] % qi
            if r != 1: good = False; break
        if good: c += 1; surv.append(S)
    cnt[k] = c
print('Part B (Fermat side only, 51 exponents): spsp(2)-type subsets by size', cnt, f'{time.time()-t0:.1f}s')
for S in surv[:40]: print('  ', S)
