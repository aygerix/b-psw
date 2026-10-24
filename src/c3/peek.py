"""Measure the D = 5 inert two-sided-smooth pool (p^2 - 1 y-smooth) by 2-adic class.
Usage: python3 peek.py y pairs_file"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from pool import build
from collections import defaultdict

y = int(sys.argv[1]); f = sys.argv[2]
rows = build(y, f)
by = defaultdict(list)
for r in rows:
    by[r['cls']].append(r)
print("y", y, "inert rows", len(rows), "max p", max(r['p'] for r in rows))
for c, R in sorted(by.items(), key=lambda t: -len(t[1])):
    qs = set()
    for r in R:
        qs |= set(r['R1']) | set(r['R2'])
    print(c, len(R), "support primes", len(qs))
