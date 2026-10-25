"""Whole-class pool metric (row primes must respect a side partition: a q dividing odd(ord) of one member and odd(z) of
another would need q | n-1 and q | n+1) for the usable DSO primes (D = 5): per 2-adic class (v2 ord, v2 z), N = #usable primes and
B = log2 phi(L1) + log2 phi(L2), L1 = lcm odd(ord), L2 = lcm odd(z) over the class (the group in which the subset
product must hit (1, -1)).  Existence needs roughly N >= B + 1; G = log2((B + 1)/N) is the cs_rank metric at full class.
Also a greedy row-pruning frontier: repeatedly drop the row prime q (on either side) whose removal costs the fewest
primes per bit saved; prints (B, N) points.  Usage: python3 dso_poolmetric.py ROWS_FILE
"""
import sys, re, math
from collections import defaultdict
from sympy import factorint, totient
rows = [tuple(map(int, t)) for t in re.findall(r'\[(\d+), (\d+), (\d+), (-?\d+)\]', open(sys.argv[1]).read())]
v2 = lambda x: (x & -x).bit_length() - 1
odd = lambda x: x >> v2(x)
cls = defaultdict(list)
for p, o, z, k in rows:
    if k == 0 or p == 2 or o % 5 == 0 or z % 5 == 0:
        continue
    if k == 1 and math.gcd(odd(o), odd(z)) != 1:
        continue
    cls[(v2(o), v2(z))].append((p, factorint(odd(o)), factorint(odd(z))))
def bits(mem):
    L1 = defaultdict(int); L2 = defaultdict(int)
    for p, f1, f2 in mem:
        for q, e in f1.items(): L1[q] = max(L1[q], e)
        for q, e in f2.items(): L2[q] = max(L2[q], e)
    b = sum(math.log2(q - 1) + (e - 1) * math.log2(q) for q, e in L1.items()) + \
        sum(math.log2(q - 1) + (e - 1) * math.log2(q) for q, e in L2.items())
    return b, L1, L2
for c in sorted(cls, key=lambda c: -len(cls[c]))[:6]:
    mem = cls[c]
    b, L1, L2 = bits(mem)
    print(f"class {c}: N = {len(mem)}, B = {b:.1f} bits (side conflicts ignored), conflicting rows {len(set(L1) & set(L2))}")
    # greedy pruning: remove the (side, q) whose primes are fewest per bit of phi(q^e)
    cur = list(mem); pts = []
    while len(cur) > 2:
        b, L1, L2 = bits(cur)
        conf = set(L1) & set(L2)
        if conf:
            # a row prime q on both sides (q | n-1 and q | n+1) is impossible: drop the side with fewer primes
            best = None
            for q in conf:
                for side in (1, 2):
                    cost = sum(1 for p, f1, f2 in cur if q in (f1 if side == 1 else f2))
                    if best is None or cost < best[0]:
                        best = (cost, side, q)
            _, side, q = best
            cur = [t for t in cur if q not in (t[1] if side == 1 else t[2])]
            continue
        pts.append((round(b, 1), len(cur)))
        best = None
        for side, L in ((1, L1), (2, L2)):
            for q, e in L.items():
                cost = sum(1 for p, f1, f2 in cur if q in (f1 if side == 1 else f2))
                gain = math.log2(q - 1) + (e - 1) * math.log2(q)
                r = cost / gain
                if best is None or r < best[0]:
                    best = (r, side, q)
        _, side, q = best
        cur = [t for t in cur if q not in (t[1] if side == 1 else t[2])]
    if not pts:
        print("   no side-consistent point"); continue
    best_ratio = max(pts, key=lambda t: t[1] / (t[0] + 1))
    sm = [t for t in pts if 60 <= t[0] <= 300]
    bs = max(sm, key=lambda t: t[1] / (t[0] + 1)) if sm else None
    print(f"   side-consistent greedy frontier: first point {pts[0]}, best N/(B+1) = {best_ratio[1]/(best_ratio[0]+1):.4f} at {best_ratio}; "
          f"best on 60 <= B <= 300: {bs} ratio {bs[1]/(bs[0]+1) if bs else 0:.4f}; points {pts[::max(1, len(pts)//8)]}")
