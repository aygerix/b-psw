"""Minimal configuration (row, side, cap) covering a list of pool primes (from a bbc 'hull primes:' line); prints the spec
string for cfg_check.py.  Usage: python3 hull_spec.py D j1 j2 OUTFILE_WITH_HULL_LINE recfile1 [recfile2 ...]"""
import sys
D, j1, j2 = map(int, sys.argv[1:4]); fn = sys.argv[4]; files = sys.argv[5:]
hull = [int(x) for l in open(fn) if l.strip().startswith('hull primes:') for x in l.split(':')[1].split()]
pool = {}
for f in files:
    for l in open(f):
        a = l.split()
        if len(a) == 7 and int(a[0]) == D and (int(a[3]), int(a[4])) == (j1, j2): pool[int(a[1])] = (int(a[5]), int(a[6]))
cap = {}
for p in hull:
    for side, x in enumerate(pool[p]):
        q = 3
        while x > 1:
            e = 0
            while x % q == 0: x //= q; e += 1
            if e:
                assert (q, 1 - side) not in cap, "row on both sides"
                cap[(q, side)] = max(cap.get((q, side), 0), e)
            q += 2
print(len(hull), ','.join(f"{q}:{'FL'[s]}:{e}" for (q, s), e in sorted(cap.items())))
