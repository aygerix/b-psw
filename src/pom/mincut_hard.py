"""Frontier of pool size against HARD bits: the part of G that is not a {2,3}-group.

G = sum over levels (q, e) of Z/phi(q^e). Its 2- and 3-primary parts admit structured subset-sum methods
(F_2 linear algebra with 2-adic lifting; the F_3 zero-sum algorithm of arXiv 2609.40321), so only the
l >= 5 part is generically hard. Hard weight of a level: e = 1: log2 of the {2,3}-free part of q - 1;
e >= 2: log2 q if q >= 5, else 0. A small multiple MU of the easy bits breaks ties.
For each lambda, the min cut maximises |S| - lambda (hard(S) + MU easy(S)). Prints the frontier with the
ranks of the 2- and 3-parts.
Usage: python3 mincut_hard.py poolfile [Tmax] [MU]
"""
import sys, math
from profile2 import primes_upto, load
from mincut import Dinic, clean_units, SCALE, INF

def split23(m):
    a = b = 0
    while m % 2 == 0: m //= 2; a += 1
    while m % 3 == 0: m //= 3; b += 1
    return a, b, m

def level_w(q, e):
    """(hard bits, easy bits) of level (q, e)."""
    if e == 1:
        a, b, r = split23(q - 1)
        return math.log2(r), a + b * math.log2(3)
    return (math.log2(q), 0.0) if q >= 5 else (0.0, math.log2(3))

def solve(rows, lam, mu):
    levels = {}
    for _, m in rows:
        for q, e in m.items():
            for ee in range(1, e + 1): levels.setdefault((q, ee), len(levels))
    n = 2 + len(rows) + len(levels); s, t = 0, 1
    D = Dinic(n)
    for i, (_, m) in enumerate(rows):
        D.add(s, 2 + i, SCALE)
        for q, e in m.items(): D.add(2 + i, 2 + len(rows) + levels[(q, e)], INF)
    for (q, e), j in levels.items():
        h, ez = level_w(q, e)
        D.add(2 + len(rows) + j, t, round(lam * SCALE * (h + mu * ez)))
        if e >= 2: D.add(2 + len(rows) + j, 2 + len(rows) + levels[(q, e - 1)], INF)
    D.maxflow(s, t)
    side = D.source_side(s)
    return [rows[i] for i in range(len(rows)) if side[2 + i]]

def describe(S):
    em = {}
    for _, m in S:
        for q, e in m.items(): em[q] = max(em.get(q, 0), e)
    hard = easy = 0.0; r2 = r3 = 0
    for q, e in em.items():
        for ee in range(1, e + 1):
            h, ez = level_w(q, ee); hard += h; easy += ez
        a, b, _ = split23(q - 1)
        r2 += a > 0; r3 += (b > 0) or (q == 3 and e >= 2)
    return hard, easy, r2, r3, len(em)

def main():
    P = primes_upto(int(sys.argv[2]) if len(sys.argv) > 2 else 2000)
    mu = float(sys.argv[3]) if len(sys.argv) > 3 else 0.01
    rows = clean_units([(p, m) for p, m in load(sys.argv[1], P)])
    print(f"pool {len(rows)} primes; MU = {mu}")
    print(" lambda |    N | hard bits | easy bits | 2-rank | 3-rank | moduli")
    seen = set()
    for lam in [0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2, 3, 5, 8, 12, 20, 40]:
        S = solve(rows, lam, mu)
        key = len(S)
        if key in seen: continue
        seen.add(key)
        h, ez, r2, r3, nm = describe(S)
        print(f" {lam:6.2f} | {len(S):4d} | {h:9.1f} | {ez:9.1f} | {r2:6d} | {r3:6d} | {nm}", flush=True)

if __name__ == "__main__":
    main()
