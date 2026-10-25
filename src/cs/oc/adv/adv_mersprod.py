#!/usr/bin/env python3
"""Adversarial re-check of Mersenne products (new code).
Fermat side: every M_q has ord 2 = q (odd), so a product n of distinct odd-q Mersenne primes is spsp(2) iff q_i | n - 1
for all i (proof: j1 = 0 for all, so need 2^odd(n-1) = 1 mod n iff q_i | odd(n-1) iff q_i | n-1).
Part A: all subsets of the 19 exponents <= 4423 by bitmask with incremental residues (Gray-code-free, plain DFS).
Part B: all subsets of size 2..5 of the 51 odd known exponents.
Part C: direct src/bpsw.py on every pair (q <= 4423) and every triple with q <= 1279 (independent of the criterion)."""
import sys, os, time, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))
from bpsw import bpsw_pass, strong_fermat2
ALL = [3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,
       44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,
       24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933,136279841]
assert len(ALL) == 51 and len(set(ALL)) == 51
t0 = time.time()
def dfs_count(Q, maxsize):
    """count subsets S (2 <= |S| <= maxsize) with prod_{j in S} M_qj == 1 mod q_i for all i in S."""
    k = len(Q)
    R = [[(pow(2, Q[j], Q[i]) - 1) % Q[i] for j in range(k)] for i in range(k)]
    hits = []; total = 0
    # residues of the running product modulo every q_i (all i, so membership test is O(|S|))
    def rec(start, S, res):
        nonlocal total
        if len(S) >= 2:
            total += 1
            if all(res[i] == 1 for i in S): hits.append(tuple(Q[i] for i in S))
        if len(S) == maxsize: return
        for j in range(start, k):
            nres = [res[i] * R[i][j] % Q[i] for i in range(k)]
            rec(j + 1, S + [j], nres)
    rec(0, [], [1] * k)
    return total, hits
QA = ALL[:19]
tot, hits = dfs_count(QA, 19)
print(f'Part A: {tot} subsets of the 19 exponents <= 4423, spsp(2) by criterion: {len(hits)} {hits[:5]}  {time.time()-t0:.1f}s', flush=True)
tot, hits = dfs_count(ALL, 5)
print(f'Part B: {tot} subsets (size 2-5) of the 51 odd exponents, spsp(2) by criterion: {len(hits)} {hits[:5]}  {time.time()-t0:.1f}s', flush=True)
# Part C: direct
npair = nb = nf = 0
for a, b in itertools.combinations(QA, 2):
    n = ((1 << a) - 1) * ((1 << b) - 1); npair += 1
    f = strong_fermat2(n); nf += f
    if f: nb += bpsw_pass(n)
print(f'Part C pairs: {npair} products, strong base-2 passes {nf}, BPSW {nb}  {time.time()-t0:.1f}s', flush=True)
Q14 = [q for q in QA if q <= 1279]
ntr = nf = 0
for S in itertools.combinations(Q14, 3):
    n = 1
    for q in S: n *= (1 << q) - 1
    ntr += 1; nf += strong_fermat2(n)
print(f'Part C triples q<=1279: {ntr} products, strong base-2 passes {nf}  {time.time()-t0:.1f}s')
