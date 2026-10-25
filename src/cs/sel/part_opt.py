"""Is the mod-4 side rule a good 'row selector'? Optimise the side partition (F, L) of the odd primes <= T0 (5 excluded)
for the usable set U(F, L) = {p : supp odd(ord_p 2) in F, supp odd(z(p)) in L}, separately for each 2-adic class.
Input: census2.c records "p j1 j2 o1 o2" (every inert p < X whose odd orders are T0-supported).
Objective: N/b of the full U (b = log2 |G| with true-order exponents); local search by single flips from the mod-4 rule
and from random starts; then the greedy frontier (frontier_sel.py) of the mod-4 rule vs the best partition.
Usage: python3 part_opt.py recfile T0 [restarts]
"""
import sys, math, random
from collections import defaultdict
sys.path.insert(0, '.')
from frontier_sel import frontier, nreq, bits
from sympy import factorint, primerange

def main():
    fn, T0 = sys.argv[1], int(sys.argv[2]); R = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    Q = [q for q in primerange(3, T0 + 1) if q != 5]
    cls = defaultdict(list)
    for l in open(fn):
        p, j1, j2, o1, o2 = map(int, l.split())
        cls[(j1, j2)].append((p, factorint(o1), factorint(o2)))
    random.seed(1)
    for c in sorted(cls, key=lambda c: -len(cls[c]))[:4]:
        rows = cls[c]
        def evalp(side):        # side[q] = 0 (F) or 1 (L)
            em = {}
            n = 0
            for p, f1, f2 in rows:
                if all(side[q] == 0 for q in f1) and all(side[q] == 1 for q in f2):
                    n += 1
                    for q, e in f1.items(): em[q] = max(em.get(q, 0), e)
                    for q, e in f2.items(): em[q] = max(em.get(q, 0), e)
            b = 1.0 + sum(bits(q, e) for q, e in em.items())
            return n, b
        def score(side):
            n, b = evalp(side); return n / b
        def local(side):
            s = score(side); improved = True
            while improved:
                improved = False
                for q in Q:
                    side[q] ^= 1; t = score(side)
                    if t > s + 1e-12: s = t; improved = True
                    else: side[q] ^= 1
            return s, dict(side)
        mod4 = {q: (0 if q % 4 == 1 else 1) for q in Q}
        n4, b4 = evalp(mod4)
        best = local(dict(mod4)); start = 'mod4'
        for r in range(R):
            s0 = {q: random.randint(0, 1) for q in Q}
            cand = local(s0)
            if cand[0] > best[0]: best = cand; start = f'random{r}'
        nb, bb = evalp(best[1])
        flips = [q for q in Q if best[1][q] != mod4[q]]
        print(f"class {c}: {len(rows)} records; mod-4 rule N = {n4}, b = {b4:.1f}, N/b = {n4 / b4:.3f}; "
              f"best partition ({start}) N = {nb}, b = {bb:.1f}, N/b = {nb / bb:.3f}; flips vs mod-4: {flips[:20]}")
        for name, side in (('mod4', mod4), ('best', best[1])):
            sub = []
            for p, f1, f2 in rows:
                if all(side[q] == 0 for q in f1) and all(side[q] == 1 for q in f2): sub.append({**f1, **f2})
            if not sub: continue
            cur = frontier(sub)
            pts = [x for x in cur if x[1] >= 2 and x[0] >= 10]
            step = max(1, len(pts) // 10)
            ex = min(((b + 1) / n, b, n) for b, n in pts) if pts else None
            print(f"   {name}: frontier", " ".join(f"({b:.0f},{n})" for b, n in pts[::step]),
                  f"| min (b+1)/N = {ex[0]:.2f} at ({ex[1]:.0f}, {ex[2]})" if ex else "")

if __name__ == "__main__":
    main()
