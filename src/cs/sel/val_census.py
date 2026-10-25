"""Independent check of census.c on a small range: recompute USABLE(T) and POOL(T) from full orders (sympy).
USABLE(T): p inert (p == +-2 mod 5), p == 3 mod 8, odd(ord_p 2) composed of primes q == 1 mod 4 (q != 5, q <= T),
odd(z(p)) composed of primes q == 3 mod 4 (q <= T).  POOL(T): (p-1)/2 and (p+1)/4 supported likewise.
Usage: python3 val_census.py X T recfile   (recfile = census.c records for [7, X) with Tmax = T)
"""
import sys
from sympy import primerange, n_order, factorint

def fibmod(k, p):
    a, b = 0, 1
    for bit in bin(k)[2:]:
        c = a * ((2 * b - a) % p) % p; d = (a * a + b * b) % p
        a, b = (d, (c + d) % p) if bit == '1' else (c, d)
    return a

def zrank(p):
    z = p + 1
    for q, e in factorint(p + 1).items():
        for _ in range(e):
            if fibmod(z // q, p) == 0: z //= q
            else: break
    return z

def odd(x):
    while x % 2 == 0: x //= 2
    return x

def main():
    X, T = int(sys.argv[1]), int(sys.argv[2])
    rec = {}
    for l in open(sys.argv[3]):
        p, o1, o2, r1, r2 = map(int, l.split()); rec[p] = (o1, o2, r1, r2)
    okF = lambda q: q % 4 == 1 and q != 5 and q <= T
    okL = lambda q: q % 4 == 3 and q <= T
    U = {}; pool = 0
    for p in primerange(7, X):
        if p % 8 != 3 or p % 5 not in (2, 3): continue
        o = odd(n_order(2, p)); z = odd(zrank(p))
        if all(okF(q) for q in factorint(o)) and all(okL(q) for q in factorint(z)):
            U[p] = (o, z)
            A, B = (p - 1) // 2, (p + 1) // 4
            if all(okF(q) for q in factorint(A)) and all(okL(q) for q in factorint(B)): pool += 1
    bad = [p for p in set(U) | set(rec) if p not in U or p not in rec or U[p] != rec[p][:2]]
    print(f"X = {X}, T = {T}: usable {len(U)} (records {len(rec)}), pool {pool}, mismatches {len(bad)} {sorted(bad)[:10]}")

if __name__ == "__main__":
    main()
