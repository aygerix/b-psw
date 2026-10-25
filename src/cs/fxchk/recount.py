import os
"""Adversarial recount of the FX pool statistics from the raw record files (new code).
Record: D p t j1 j2 o1 o2 (t = 0 inert, 1 split)."""
import sys
from collections import defaultdict, Counter
B = os.path.expanduser('~/b-psw/src/cs/')
F200 = [B + 'e1/c3_1e10.txt', B + 'e1/chm200_pool.txt', B + 'fx/c3_xD_1e10.txt', B + 'fx/chm200_pool_xD.txt']
def load(files):
    pool = {}; conflict = 0; dup = 0
    for fn in files:
        with open(fn) as f:
            for l in f:
                a = l.split()
                if len(a) != 7: continue
                r = tuple(int(x) for x in a)
                k = (r[0], r[1])
                if k in pool:
                    dup += 1
                    if pool[k] != r: conflict += 1
                pool[k] = r
    return pool, dup, conflict
pool, dup, conf = load(F200)
print(f"T0=200 union: {len(pool)} (D,p) records, duplicates {dup}, conflicting duplicates {conf}")
Ds = sorted({k[0] for k in pool}, key=lambda d: (abs(d), d))
cls = defaultdict(Counter); comp = defaultdict(Counter); ij2 = defaultdict(Counter)
for (D, p), r in pool.items():
    _, _, t, j1, j2, o1, o2 = r
    cls[D][(j1, j2)] += 1; comp[(D, j1, j2)][t] += 1
    if t == 0 and j1 == 2: ij2[D][min(j2, 2)] += 1
print("D: records, #classes, inert j1=2 (j2=0/1/>=2), total inert j1=2, class (2,0) size [inert, split], largest class")
allc = []
for D in Ds:
    c = cls[D]; big = c.most_common(1)[0]
    for k, v in c.items(): allc.append((v, D, k))
    print(f"  {D:4d}: {sum(c.values()):6d} {len(c):3d}  {ij2[D][0]:5d}/{ij2[D][1]:5d}/{ij2[D][2]:3d}  tot {ij2[D][0]+ij2[D][1]+ij2[D][2]:5d}  (2,0) {c[(2,0)]:5d} [{comp[(D,2,0)][0]},{comp[(D,2,0)][1]}]  largest {big}")
allc.sort(reverse=True); print("top classes over the 24 D:", allc[:6])
print("number of (D,class) instances:", len(allc), " ten D:", sum(1 for v, D, k in allc if D in (5,-7,-11,13,-15,17,-19,-23,29,-31)))
# 5-adic cells, class (2,0) inert, D = 5 and -15
def leg(a, p):
    r = pow(a % p, (p - 1) // 2, p); return 0 if r == 0 else (1 if r == 1 else -1)
for D in (5, -15):
    cell = Counter()
    for (DD, p), r in pool.items():
        if DD == D and r[2] == 0 and r[3] == 2 and r[4] == 0: cell[(p % 3, leg(p, 5))] += 1
    print(f"5-adic cells D={D} class (2,0) inert:", dict(cell))
# per-bin pools
bins = [0, 1e7, 1e8, 1e9, 1e10, 1e300]
def binned(D, c):
    h = [0] * (len(bins) - 1)
    for (DD, p), r in pool.items():
        if DD == D and (r[3], r[4]) == c:
            for i in range(len(bins) - 1):
                if bins[i] <= p < bins[i + 1]: h[i] += 1
    return h
h15 = binned(-15, (2, 0)); h11 = binned(-11, (2, 0)); h5 = binned(5, (2, 0))
print("bins", bins[1:-1], "\n (2,0) counts D=-15", h15, " D=-11", h11, " D=5", h5)
print(" ratio -15/-11:", [round(a / b, 2) if b else None for a, b in zip(h15, h11)], " ratio -15/5:", [round(a / b, 2) if b else None for a, b in zip(h15, h5)])
# T0 = 300 superset check for D = -15 (2,0) etc.
p300, d3, c3 = load([B + 'fx/c3_T300_1e10.txt'])
for D, c in ((-15, (2, 0)), (-7, (1, 2)), (5, (2, 0))):
    s200 = {p for (DD, p), r in pool.items() if DD == D and (r[3], r[4]) == c}
    s300 = {p for (DD, p), r in p300.items() if DD == D and (r[3], r[4]) == c}
    miss = s200 - s300
    print(f"T0=300 census-only pool D={D} {c}: {len(s300)}; T0=200 pool {len(s200)}; T0=200 primes missing from T0=300 census: {len(miss)}"
          f" (of which p > 1e10: {sum(1 for p in miss if p > 1e10)}); max p in T0=200 pool {max(s200)}")
