"""Pairwise-closed graph on inert usable primes inside each 2-adic class (audit rows "p o z s"):
p ~ q iff odd(o_q) | p-1, odd(o_p) | q-1, odd(z_q) | p+1, odd(z_p) | q+1.  usage: python3 aud_closed.py D rows..."""
import sys, math
from collections import defaultdict
v2 = lambda x: (x & -x).bit_length() - 1
odd = lambda x: x >> v2(x)
D = int(sys.argv[1])
for fn in sys.argv[2:]:
    R = [tuple(map(int, l.split())) for l in open(fn)]
    cls = defaultdict(list)
    for p, o, z, s in R:
        if s != -1: continue
        if D == 5 and (o % 5 == 0 or z % 5 == 0): continue
        cls[(v2(o), v2(z))].append((p, odd(o), odd(z)))
    edges = 0; pairs = 0
    for c, m in cls.items():
        for i in range(len(m)):
            p, a, b = m[i]
            for j in range(i):
                q, c2, d2 = m[j]; pairs += 1
                if (p - 1) % c2 == 0 and (q - 1) % a == 0 and (p + 1) % d2 == 0 and (q + 1) % b == 0: edges += 1
    print(f"{fn}: inert usable {sum(map(len, cls.values()))}, same-class pairs {pairs}, closed edges {edges}")
