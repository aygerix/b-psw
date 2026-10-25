"""FX (1): decompose the E1 smooth pool (census3 p < 1e10 + CHM y = 200) by 2-adic class and by small-prime local data.
Usage: python3 pooldec.py rec1 [rec2 ...]
Checks R1 on every inert j1 = 2 record (j2 = 0 iff (Q/p) = +1), then prints:
  (A) per D: inert j1 = 2 records split j2 = 0 / 1  (the (2,0) / (2,1) pools minus split primes)
  (B) per-bin pool of (2,0) for D = 5, -15, 17 and of (2,0) + (2,1) for the other D
  (C) D = 5 versus D = -15, class (2,0), split by (p mod 3, (p/5)): the 5-adic test."""
import sys
from collections import defaultdict
def leg(a, p):
    r = pow(a % p, (p - 1) // 2, p)
    return 0 if r == 0 else (1 if r == 1 else -1)
recs = []
for fn in sys.argv[1:]:
    for l in open(fn):
        a = l.split()
        if len(a) == 7: recs.append(tuple(int(x) for x in a))
# dedup (D, p): census and CHM overlap below 1e10
seen = {}
for r in recs: seen[(r[0], r[1])] = r
recs = list(seen.values())
print("records (dedup):", len(recs))
DL = sorted({r[0] for r in recs}, key=lambda d: (abs(d), d))
viol = 0; A = defaultdict(lambda: [0, 0, 0])
for D, p, t, j1, j2, o1, o2 in recs:
    if t == 0 and j1 == 2:
        Q = (1 - D) // 4
        if (j2 == 0) != (leg(Q, p) == 1): viol += 1
        A[D][min(j2, 2)] += 1
print("R1 violations on inert j1 = 2 records:", viol)
print("(A) inert j1 = 2 records: j2 = 0 / 1 / >= 2, share(j2 = 0), total")
for D in DL:
    a = A[D]; print(f"   D = {D:4d} Q = {(1 - D) // 4:3d}: {a[0]:5d} {a[1]:5d} {a[2]:3d}  {a[0] / max(1, sum(a)):.3f}  {sum(a)}")
bins = [(0, 1e7), (1e7, 1e8), (1e8, 1e9), (1e9, 1e10), (1e10, 1e40)]
print("(B) per p-bin: inert+split class (2,0) [and (2,0)+(2,1)] pool sizes")
B = defaultdict(lambda: [0] * len(bins)); B2 = defaultdict(lambda: [0] * len(bins))
for D, p, t, j1, j2, o1, o2 in recs:
    k = next(i for i, (lo, hi) in enumerate(bins) if lo <= p < hi)
    if (j1, j2) == (2, 0): B[D][k] += 1
    if j1 == 2 and j2 in (0, 1): B2[D][k] += 1
print("   bins:", ["%.0e-%.0e" % b for b in bins])
for D in DL:
    print(f"   D = {D:4d}: (2,0) {B[D]}   (2,0)+(2,1) {B2[D]}   ratio -15/(this) (2,0): " +
          " ".join("%.2f" % (B[-15][k] / B[D][k]) if B[D][k] else "-" for k in range(len(bins))))
print("(C) class (2,0) inert, split by (p mod 3, (p/5)), D = 5 vs D = -15 (and D = 17, -7 j1 = 2 inert for reference)")
for D in (5, -15, 17, -7, -11, -19, -23, -31, 13, 29):
    c = defaultdict(int)
    for D2, p, t, j1, j2, o1, o2 in recs:
        if D2 == D and t == 0 and j1 == 2 and (D in (5, -15, 17) and j2 == 0 or D not in (5, -15, 17) and j2 in (0, 1)):
            c[(p % 3, leg(5, p) if p != 5 else 0)] += 1
    print(f"   D = {D:4d}: " + "  ".join(f"p%3={k[0]},(p/5)={k[1]:+d}: {c[k]:5d}" for k in sorted(c)))
