"""Build Pomerance-type pools from smooth-neighbor pairs and optimise the A/B partition.
For x with x, x+1 z-smooth: p = 2x+1. Keep p prime, 5 not | x(x+1)  (=> (5/p) = -1, inert).
Requirement sets (exact, via Lemma 1 using actual orders):
   R1(p) = prime-power factorisation of odd(ord_p 2)   (must divide L1, n = 1 mod L1)
   R2(p) = prime-power factorisation of odd(z(p))      (must divide L2, n = -1 mod L2)
Class = (v2(ord_p 2), v2(z(p))).  Objective: N - log2|G|,  G = (Z/L1)^* x (Z/L2)^* x Z/2."""
import sys, math, random, json
from collections import defaultdict
from sympy import isprime, factorint, primerange
import numpy as np

def v2(x): return (x & -x).bit_length() - 1
def odd(x): return x >> v2(x)

def order_mod(a, p, fac):           # fac = factorint(p-1)
    o = p - 1
    for q, e in fac.items():
        for _ in range(e):
            if pow(a, o // q, p) == 1: o //= q
            else: break
    return o

def fib_mod(k, p):
    a, b = 0, 1
    for bit in bin(k)[2:]:
        c = a * ((2 * b - a) % p) % p; d = (a * a + b * b) % p
        a, b = (d, (c + d) % p) if bit == '1' else (c, d)
    return a

def zentry(p, fac):                 # fac = factorint(p+1) (inert p)
    o = p + 1
    for q, e in fac.items():
        for _ in range(e):
            if fib_mod(o // q, p) == 0: o //= q
            else: break
    return o

def build(z, fname):
    xs = [int(l) for l in open(fname)]
    rows = []
    for x in xs:
        p = 2 * x + 1
        if x % 5 == 0 or (x + 1) % 5 == 0 or p == 5 or not isprime(p): continue
        fm, fp = factorint(p - 1), factorint(p + 1)
        o2 = order_mod(2, p, fm); zz = zentry(p, fp)
        rows.append(dict(p=p, cls=(v2(o2), v2(zz)), R1=factorint(odd(o2)), R2=factorint(odd(zz)),
                         F1=factorint(odd(p - 1)), F2=factorint(odd(p + 1))))
    return rows

def evaluate(rows, assign, use_orders=True):
    """assign: dict q -> 'A'/'B'/'X'. Returns N, log2|G|, L1 exps, L2 exps, pool."""
    pool, e1, e2 = [], defaultdict(int), defaultdict(int)
    for r in rows:
        R1, R2 = (r['R1'], r['R2']) if use_orders else (r['F1'], r['F2'])
        if r['p'] in assign and assign[r['p']] != 'X': continue
        if all(assign.get(q) == 'A' for q in R1) and all(assign.get(q) == 'B' for q in R2):
            pool.append(r)
            for q, e in R1.items(): e1[q] = max(e1[q], e)
            for q, e in R2.items(): e2[q] = max(e2[q], e)
    lg = 1.0 + sum(math.log2(q - 1) + (e - 1) * math.log2(q) for q, e in list(e1.items()) + list(e2.items()))
    return len(pool), lg, dict(e1), dict(e2), pool

def anneal(rows, qs, iters=20000, T0=3.0, seed=0, use_orders=True):
    rnd = random.Random(seed)
    assign = {q: rnd.choice('AB') for q in qs}
    N, lg, *_ = evaluate(rows, assign, use_orders); cur = N - lg
    best = (cur, dict(assign), N, lg)
    for it in range(iters):
        T = T0 * (1 - it / iters) + 0.05
        q = rnd.choice(qs); old = assign[q]
        assign[q] = rnd.choice([c for c in 'ABX' if c != old])
        N, lg, *_ = evaluate(rows, assign, use_orders); val = N - lg
        if val >= cur or rnd.random() < math.exp((val - cur) / T):
            cur = val
            if val > best[0]: best = (val, dict(assign), N, lg)
        else:
            assign[q] = old
    return best

if __name__ == "__main__":
    z = int(sys.argv[1]); fname = sys.argv[2]
    rows = build(z, fname)
    bycls = defaultdict(list)
    for r in rows: bycls[r['cls']].append(r)
    print(f"z={z}: inert primes p=2x+1 from pairs: {len(rows)}; max p ~ 2^{max(r['p'] for r in rows).bit_length()}")
    for c in sorted(bycls, key=lambda c: -len(bycls[c]))[:8]:
        print("  class (j1,j2) =", c, "size", len(bycls[c]))
    qs = [q for q in primerange(3, z + 1) if q != 5]
    out = {}
    for c in sorted(bycls, key=lambda c: -len(bycls[c]))[:5]:
        R = bycls[c]
        for uo in (False, True):
            best = max((anneal(R, qs, iters=6000, seed=s, use_orders=uo) for s in range(4)), key=lambda b: b[0])
            print(f"class {c} orders={uo}: best N - log2|G| = {best[0]:.1f}  (N={best[2]}, log2|G|={best[3]:.1f})"
                  f"  A={[q for q in qs if best[1][q]=='A']} B={[q for q in qs if best[1][q]=='B']}")
            out[f"{c}-{uo}"] = dict(val=best[0], N=best[2], lg=best[3], assign=best[1])
    json.dump(out, open(f"../data/pool_best_z{z}.json", "w"), default=str, indent=1)
