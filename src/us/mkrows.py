"""Convert usplit3 output ('D p ord z' per usable split prime) into a GP ROWS file, keeping only
admissible primes (Lemma 17a: gcd(odd ord, odd z) = 1 and min(v2 ord, v2 z) <= 1).
Usage: python3 mkrows.py in.txt out.gp"""
import sys
from math import gcd
def v2(x): return (x & -x).bit_length() - 1
def odd(x): return x >> v2(x)
rows = []
for l in open(sys.argv[1]):
    f = l.split()
    if len(f) != 4: continue
    D, p, o, z = map(int, f)
    if gcd(odd(o), odd(z)) == 1 and min(v2(o), v2(z)) <= 1:
        rows.append((D, p, o, z))
open(sys.argv[2], 'w').write("ROWS=[" + ",".join(f"[{D},{p},{o},{z}]" for D, p, o, z in rows) + "];\n")
from collections import Counter
print(len(rows), "admissible rows;", dict(Counter(r[0] for r in rows)))
for D in sorted(set(r[0] for r in rows)):
    print(D, [r[1] for r in rows if r[0] == D])
