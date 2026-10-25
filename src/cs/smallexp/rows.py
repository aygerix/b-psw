#!/usr/bin/env python3
"""Designed row sets. pierpont: q - 1 = 2^a 3^b; hamming: q - 1 = 2^a 3^b 5^c. 5 excluded (D = 5 inert pools).
Usage: python3 rows.py pierpont|hamming Y   -> space-separated primes (2 included)"""
import sys
from sympy import isprime
def smooth_vals(Y, ps):
    out = {1}
    for p in ps:
        new = set()
        for v in out:
            w = v
            while w <= Y:
                new.add(w); w *= p
        out = new
    return out
def rows(kind, Y):
    ps = [2, 3] if kind == 'pierpont' else [2, 3, 5]
    return sorted({2} | {d + 1 for d in smooth_vals(Y, ps) if d + 1 != 5 and isprime(d + 1)})
if __name__ == "__main__":
    print(" ".join(map(str, rows(sys.argv[1], int(sys.argv[2])))))
