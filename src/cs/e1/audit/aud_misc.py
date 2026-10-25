"""Misc audit numbers: per-prime own-requirement weight w_p = sum over q^a || o1 o2 of log2 phi(q^a) (median per D, all
classes, and admissible classes only); B_all = total item cost of the D = 5 (2,0) instance and the LP value |P| b / B_all.
Usage: python3 aud_misc.py c3file chmfile inst_D5_2_0"""
import sys, math, statistics
from collections import defaultdict
PR = [q for q in range(3, 201) if all(q % d for d in range(2, int(q ** .5) + 1))]
def fac(x):
    f = {}
    for q in PR:
        while x % q == 0: x //= q; f[q] = f.get(q, 0) + 1
    return f
rec = {}
for fn in sys.argv[1:3]:
    for l in open(fn):
        t = tuple(map(int, l.split()))
        if len(t) == 7: rec[(t[0], t[1])] = t
w = defaultdict(list); wa = defaultdict(list)
for (D, p), t in rec.items():
    s = sum(math.log2(q - 1) + (a - 1) * math.log2(q) for o in (t[5], t[6]) for q, a in fac(o).items())
    w[D].append(s)
    if min(t[3], t[4]) <= 1: wa[D].append(s)
print("median w_p per D (all / admissible):", {D: (round(statistics.median(w[D]), 1), round(statistics.median(wa[D]), 1)) for D in sorted(w)})
f = open(sys.argv[3]); nP, nI, nR = map(int, f.readline().split())
B = sum(float(f.readline().split()[4]) for _ in range(nI))
print(f"D=5 (2,0): |P| = {nP}, B_all = {B:.3f}; |P| b / B_all at b = 60, 140, 214, 275:",
      [round(nP * b / B, 2) for b in (60, 140, 214, 275)], "; |P|/2 =", nP / 2)
# straggler fraction per decade (all ten D pooled): census3 records whose (D, p) is absent from the CHM records
c3 = set(); chm = set()
for fn, dst in ((sys.argv[1], c3), (sys.argv[2], chm)):
    for l in open(fn):
        a = l.split(); dst.add((int(a[0]), int(a[1])))
dc, ds, dh = defaultdict(int), defaultdict(int), defaultdict(int)
for (D, p) in c3:
    k = len(str(p)) - 1; dc[k] += 1
    if (D, p) not in chm: ds[k] += 1
for (D, p) in chm: dh[len(str(p)) - 1] += 1
print("decade: stragglers / CHM (fraction):", "  ".join(f"1e{k}: {ds[k]}/{dh[k]} ({ds[k] / dh[k]:.3%})" for k in range(4, 10)))
s_const = (ds[8] + ds[9]) / (dh[8] + dh[9])
# linear trend of the fraction in the decade index over 1e6..1e9, extrapolated per decade above 1e10
xs = [6, 7, 8, 9]; ys = [ds[k] / dh[k] for k in xs]
mx, my = sum(xs) / 4, sum(ys) / 4
sl = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
above = {k: v for k, v in dh.items() if k >= 10}
est_const = s_const * sum(above.values())
est_trend = sum(v * (my + sl * (k - mx)) for k, v in above.items())
print(f"CHM records above 1e10 (all D): {sum(above.values())}; extrapolated missing stragglers: constant rate {s_const:.4f} -> "
      f"{est_const:.0f} (line's method), linear trend (slope {sl:.4f}/decade) -> {est_trend:.0f}")
