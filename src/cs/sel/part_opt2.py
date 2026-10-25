"""Small-b existence gap, optimised directly: for each 2-adic class, local search over side partitions (F, L) of the odd
primes <= T0 maximising  f = max over greedy-frontier points with b <= BMAX of N/(b+1)  (N >= b+1 is needed for a
subset-product solution to exist; b <= ~250 is the range where density-1 algorithms cost <= ~2^73).
Input: census2.c records "p j1 j2 o1 o2".  Usage: python3 part_opt2.py recfile T0 BMAX restarts [BMIN]
"""
import sys, random
from collections import defaultdict
sys.path.insert(0, '.')
from frontier_sel import frontier
from sympy import factorint, primerange

def main():
    fn, T0, BMAX, R = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]); BMIN = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
    Q = [q for q in primerange(3, T0 + 1) if q != 5]
    cls = defaultdict(list)
    for l in open(fn):
        p, j1, j2, o1, o2 = map(int, l.split())
        cls[(j1, j2)].append((p, factorint(o1), factorint(o2)))
    random.seed(7)
    for c in sorted(cls, key=lambda c: -len(cls[c]))[:3]:
        rows = cls[c]
        def f(side):
            sub = [{**f1, **f2} for p, f1, f2 in rows
                   if all(side[q] == 0 for q in f1) and all(side[q] == 1 for q in f2)]
            if not sub: return 0.0, None
            best = (0.0, None)
            for b, n in frontier(sub):
                if BMIN <= b <= BMAX and n / (b + 1) > best[0]: best = (n / (b + 1), (round(b, 1), n))
            return best
        def local(side):
            s = f(side); improved = True
            while improved:
                improved = False
                for q in random.sample(Q, len(Q)):
                    side[q] ^= 1; t = f(side)
                    if t[0] > s[0] + 1e-12: s = t; improved = True
                    else: side[q] ^= 1
            return s, dict(side)
        mod4 = {q: (0 if q % 4 == 1 else 1) for q in Q}
        res = [('mod4', f(mod4))]
        best = local(dict(mod4)); tag = 'from mod4'
        for r in range(R):
            cand = local({q: random.randint(0, 1) for q in Q})
            if cand[0][0] > best[0][0]: best = cand; tag = f'random{r}'
        print(f"class {c} ({len(rows)} records), {BMIN} <= b <= {BMAX}: mod-4 best N/(b+1) = {res[0][1][0]:.4f} at {res[0][1][1]}; "
              f"optimised ({tag}) = {best[0][0]:.4f} at {best[0][1]}  -> existence gap {1 / best[0][0]:.2f}", flush=True)

if __name__ == "__main__":
    main()
