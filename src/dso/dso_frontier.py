# Frontier (log2|G|, N) for the doubly-small-order pool (ord_p 2, z(p) <= 1500), all classes.
import sys, re, math, random
from collections import defaultdict
from sympy import factorint
sys.path.insert(0, '.')
from pool import evaluate
from frontier import anneal_l
txt = open('dso/dso_primes.txt').read()
L = [tuple(map(int, t)) for t in re.findall(r'\[(\d+), (\d+), (\d+), (-?\d+)\]', txt)]
v2 = lambda x: (x & -x).bit_length() - 1
odd = lambda x: x >> v2(x)
rows = [dict(p=p, cls=(v2(o), v2(z)), R1=factorint(odd(o)), R2=factorint(odd(z)), jac=j) for p, o, z, j in L if p != 5]
bycls = defaultdict(list)
for r in rows: bycls[r['cls']].append(r)
print("classes:", sorted(((c, len(v)) for c, v in bycls.items()), key=lambda t: -t[1])[:8])
for c in sorted(bycls, key=lambda c: -len(bycls[c]))[:4]:
    qs = sorted({q for r in bycls[c] for q in list(r['R1']) + list(r['R2'])})
    pts = []
    for lam in (0.0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75):
        b = max((anneal_l(bycls[c], qs, lam, 6000, s) for s in range(3)), key=lambda t: t[0])
        pts.append((round(b[3], 1), b[2]))
    print(f"class {c} |class|={len(bycls[c])} #q={len(qs)} frontier (log2|G|, N):", sorted(set(pts)))
