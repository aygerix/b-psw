"""Cross-checks of the E2 census against Attempt D (src/dso/dso_primes.txt) and census2 records (c2_1e10_200.txt),
and against the direct sieve output (direct_1e9_X6e4.txt) filtered to smaller X.  Uses the audit's own recomputed rows
(rows_X*.txt, "p o z s").  usage: python3 aud_xcheck.py"""
import re, math
def rows(fn):
    return [tuple(map(int, l.split())) for l in open(fn) if l.strip()]
R = {X: rows(f"rows_X{X}.txt") for X in (10000, 20000, 60000)}
v2 = lambda x: (x & -x).bit_length() - 1
odd = lambda x: x >> v2(x)
# Attempt D
att = [tuple(map(int, t)) for t in re.findall(r'\[(\d+), (\d+), (\d+), (-?\d+)\]', open('../../../dso/dso_primes.txt').read())]
attS = {r[0]: r for r in att}
c1500 = [r for r in R[10000] if r[1] <= 1500 and r[2] <= 1500]
print(f"Attempt D file: {len(att)} rows (min p {min(attS)}, max p {max(attS)}); census rows with ord, z <= 1500: {len(c1500)}")
cS = {r[0]: r for r in c1500}
print("  in census not in Attempt D:", sorted(set(cS) - set(attS)), " in Attempt D not in census:", sorted(set(attS) - set(cS)))
print("  rows differing:", [p for p in set(cS) & set(attS) if cS[p] != attS[p]][:10])
# direct sieve filtered
dr = rows('../direct_1e9_X6e4.txt')
for X in (10000, 20000, 60000):
    d = {(p, o, z, s) for p, o, z, s in dr if o <= X and z <= X}
    c = set(R[X])
    print(f"direct (p <= 1e9, ord, z <= 6e4) filtered to X = {X}: {len(d)} rows; census {len(c)}; identical {d == c}")
# census2: inert, 7 <= p < 1e10, odd(ord), odd(z) supported on odd primes <= 200 other than 5
c2 = {}
for l in open('c2_1e10_200.txt'):
    p, j1, j2, o1, o2 = map(int, l.split()); c2[p] = (j1, j2, o1, o2)
def smooth200(m):
    for q in range(3, 201, 2):
        if q == 5: continue
        while m % q == 0: m //= q
    return m == 1
for X in (10000, 20000, 60000):
    sel = [r for r in R[X] if r[3] == -1 and smooth200(odd(r[1])) and smooth200(odd(r[2]))]
    inC2 = [r for r in sel if r[0] in c2]
    same = [r for r in inC2 if c2[r[0]] == (v2(r[1]), v2(r[2]), odd(r[1]), odd(r[2]))]
    print(f"X = {X}: inert DSO primes with 200-smooth odd orders (5 excluded): {len(sel)}; in census2 {len(inC2)}; identical orders {len(same)}; "
          f"missing {[r for r in sel if r[0] not in c2]}")
