"""List the pool primes usable in a configuration and write a PARI script that re-verifies them independently with
e1/audit/verify_cfg.gp (orders recomputed from scratch, bits = log2 phi(E1 E2)).
Usage: python3 cfg_check.py D j1 j2 "q:s:e,..." OUT.gp recfile1 [recfile2 ...]"""
import sys, math
D, j1, j2 = map(int, sys.argv[1:4]); spec = sys.argv[4]; out = sys.argv[5]; files = sys.argv[6:]
cap = {0: {}, 1: {}}
for t in spec.split(','):
    q, s, e = t.split(':'); cap[0 if s == 'F' else 1][int(q)] = int(e)
def ok(x, side):
    for q in list(range(3, 400)):
        if x % q == 0 and all(q % d for d in range(2, int(q ** .5) + 1)):
            e = 0
            while x % q == 0: x //= q; e += 1
            if cap[side].get(q, 0) < e: return False
    return x == 1
pool = {}
for fn in files:
    for l in open(fn):
        a = l.split()
        if len(a) == 7 and int(a[0]) == D and (int(a[3]), int(a[4])) == (j1, j2): pool[int(a[1])] = (int(a[5]), int(a[6]))
use = sorted(p for p, (o1, o2) in pool.items() if ok(o1, 0) and ok(o2, 1))
bits = sum(math.log2(q - 1) + (e - 1) * math.log2(q) for s in (0, 1) for q, e in cap[s].items())
print(f"D = {D} class ({j1},{j2}): pool {len(pool)}, usable {len(use)}, bits {bits:.3f}, N - 0.25 bits = {len(use) - 0.25 * bits:.4f}, "
      f"N - 0.2 bits = {len(use) - 0.2 * bits:.4f}")
E1 = [[q, e] for q, e in sorted(cap[0].items())]; E2 = [[q, e] for q, e in sorted(cap[1].items())]
with open(out, 'w') as f:
    f.write('read("../e1/audit/verify_cfg.gp");\n')
    f.write(f"vcfg({D}, {j1}, {j2}, {E1}, {E2}, {use});\n")
