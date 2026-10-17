"""Variant: Selfridge D=-7 (requires (5/n)=+1), P=1, Q=2. Inert p: (p/7)=-1. 5 is now usable on either side."""
import sys, math, random, json
from collections import defaultdict
from sympy import isprime, factorint, primerange, n_order
from bpsw import lucas_UV, jacobi, strong_fermat2, strong_lucas
import pool as P0

def v2(x): return (x & -x).bit_length() - 1
def odd(x): return x >> v2(x)

def zq2(p):
    o = p + 1
    for q, e in factorint(p + 1).items():
        for _ in range(e):
            if lucas_UV(o // q, 1, 2, p)[0] == 0: o //= q
            else: break
    return o

def build7(fname):
    rows = []
    for x in (int(l) for l in open(fname)):
        p = 2 * x + 1
        if p in (7,) or p % 7 not in (3, 5, 6) or not isprime(p): continue
        o2 = n_order(2, p); zz = zq2(p)
        rows.append(dict(p=p, cls=(v2(o2), v2(zz)), R1=factorint(odd(o2)), R2=factorint(odd(zz)),
                         F1=factorint(odd(p - 1)), F2=factorint(odd(p + 1))))
    return rows

if __name__ == "__main__":
    # sanity: criterion vs direct test for D=-7 composites
    random.seed(2); Ps = [p for p in primerange(11, 3000)]
    ag = dis = 0
    for _ in range(30000):
        fs = random.sample(Ps, random.choice([2, 3, 3, 5])); n = math.prod(fs)
        if jacobi(5, n) != 1 or jacobi(-7, n) != -1: continue
        s, t = v2(n - 1), v2(n + 1)
        j1 = {v2(n_order(2, p)) for p in fs}; j2 = {v2(zq2(p)) if p % 7 in (3,5,6) else None for p in fs}
        pred = (None not in j2 and len(j1) == 1 and len(j2) == 1 and max(j1) <= s and max(j2) <= t and
                all((n - 1) % odd(n_order(2, p)) == 0 and (n + 1) % odd(zq2(p)) == 0 for p in fs))
        # (split primes w.r.t. -7 handled by direct test only; pred False for them)
        d = strong_fermat2(n) and strong_lucas(n)
        ag += (pred == d); dis += (pred != d)
    print("D=-7 criterion check: agree", ag, "disagree", dis)
    z = int(sys.argv[1]); rows = build7(sys.argv[2])
    bycls = defaultdict(list)
    for r in rows: bycls[r['cls']].append(r)
    print(f"z={z}: D=-7 inert primes from pairs: {len(rows)};", {c: len(v) for c, v in sorted(bycls.items(), key=lambda kv: -len(kv[1]))[:6]})
    from frontier import anneal_l
    qs = [q for q in primerange(3, z + 1) if q != 7]
    for c in sorted(bycls, key=lambda c: -len(bycls[c]))[:3]:
        pts = set()
        for lam in (0.05, 0.1, 0.2, 0.35, 0.5):
            b = max((anneal_l(bycls[c], qs, lam, int(sys.argv[3]), s) for s in range(3)), key=lambda t: t[0])
            pts.add((round(b[3] + 1, 1), b[2], round(b[2] - b[3] - 1, 1)))   # +1 bit for (5/n)=+1 parity
        print(f"z={z} D=-7 class {c} frontier (log2|G|, N, N-log2|G|):", sorted(pts))
