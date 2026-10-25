"""Row-prime multiplicity in the usable D = 5 DSO pool: for each odd row prime q (q | odd(ord) -> Fermat side,
q | odd(z) -> Lucas side), the number of usable primes requiring it; summary for q <= 200 and q > 200.
usage: python3 aud_rows.py rows_X20000.txt"""
import sys, math
from collections import Counter, defaultdict
from sympy import primefactors
v2 = lambda x: (x & -x).bit_length() - 1
odd = lambda x: x >> v2(x)
R = [tuple(map(int, l.split())) for l in open(sys.argv[1])]
U = [r for r in R if r[3] != 0 and r[1] % 5 and r[2] % 5 and (r[3] == -1 or math.gcd(odd(r[1]), odd(r[2])) == 1)]
cnt = Counter()
for p, o, z, s in U:
    for q in primefactors(odd(o)): cnt[('F', q)] += 1
    for q in primefactors(odd(z)): cnt[('L', q)] += 1
small = {k: v for k, v in cnt.items() if k[1] <= 200}; big = {k: v for k, v in cnt.items() if k[1] > 200}
mb = Counter(big.values())
print(f"usable {len(U)}; (side,row) pairs: q <= 200: {len(small)} (mean multiplicity {sum(small.values())/len(small):.2f}); "
      f"q > 200: {len(big)} (mean multiplicity {sum(big.values())/len(big):.3f}); multiplicity histogram q > 200: {sorted(mb.items())[:8]}")
need = sum(1 for p, o, z, s in U if all(q <= 200 for q in primefactors(odd(o) * odd(z))))
print(f"usable primes with all rows <= 200: {need}")
