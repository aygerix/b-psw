#!/usr/bin/env python3
"""Control for the Lemma-1 prefilter used in oc_post.subset_search: split it into its Fermat and Lucas halves and
compare with src/bpsw.py strong_fermat2 / strong_lucas on all products of 2 distinct primes < 3000 and random triples.
Half-criteria (per prime, with the Selfridge D of n): one 2-adic class j1 (resp. j2), odd orders | n - 1 (resp. n + 1),
v2 room."""
import sys, os, itertools, random, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from bpsw import strong_fermat2, strong_lucas, selfridge_D, lucas_UV
from sympy import primerange, n_order, factorint
P = list(primerange(3, 3000))
def v2(x): return (x & -x).bit_length() - 1
def zrank(p, D):
    Q = (1 - D) // 4; eps = pow(D % p, (p - 1)//2, p); eps = -1 if eps == p - 1 else eps
    N = p - eps; z = N
    for r, e in factorint(N).items():
        for _ in range(e):
            if z % r == 0 and lucas_UV(z // r, 1, Q, p)[0] == 0: z //= r
            else: break
    return z
def pred(ps):
    n = math.prod(ps); D = selfridge_D(n)
    if D is None: return None, None
    Q = (1 - D)//4
    if any((D * Q * 2) % p == 0 for p in ps): return None, None
    s, t = v2(n - 1), v2(n + 1)
    o = [n_order(2, p) for p in ps]; z = [zrank(p, D) for p in ps]
    jf = set(v2(x) for x in o); jl = set(v2(x) for x in z)
    F = len(jf) == 1 and all((n - 1) % (x >> v2(x)) == 0 for x in o) and min(jf) <= s
    L = len(jl) == 1 and all((n + 1) % (x >> v2(x)) == 0 for x in z) and min(jl) <= t
    return F, L
cnt = dict(n=0, fa=0, la=0, fh=0, lh=0)
def run(ps):
    F, L = pred(ps)
    if F is None: return
    n = math.prod(ps); a, b = strong_fermat2(n), strong_lucas(n)
    cnt['n'] += 1; cnt['fa'] += (F == a); cnt['la'] += (L == b); cnt['fh'] += a; cnt['lh'] += b
    if F != a or L != b: print('DISAGREE', ps, F, a, L, b)
for a, b in itertools.combinations(P, 2): run((a, b))
random.seed(7)
for _ in range(100000): run(tuple(random.sample(P[:200], 3)))
print(f"tested {cnt['n']}: Fermat-half agree {cnt['fa']}, Lucas-half agree {cnt['la']}; true spsp(2) {cnt['fh']}, true strong Lucas psp {cnt['lh']}")
