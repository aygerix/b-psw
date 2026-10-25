"""Independent brute-force validator for census3.c on small X: exact ord_p 2 (sympy n_order) and exact rank of
apparition z_D(p) (smallest divisor d of p - (D/p) with U_d == 0, U from 2x2 matrix powers), then the same filter.
Usage: python3 val3.py X T0 census3_recfile   -> prints mismatches and counts per D."""
import sys
from sympy import primerange, n_order, factorint, jacobi_symbol, divisors

DL = [5, -7, -11, 13, -15, 17, -19, -23, 29, -31]

def matU(k, Q, p):           # U_k(1, Q) mod p via [[1, -Q], [1, 0]]^k ; U_k = M^k[1][0]
    def mul(A, B):
        return [[(A[0][0] * B[0][0] + A[0][1] * B[1][0]) % p, (A[0][0] * B[0][1] + A[0][1] * B[1][1]) % p],
                [(A[1][0] * B[0][0] + A[1][1] * B[1][0]) % p, (A[1][0] * B[0][1] + A[1][1] * B[1][1]) % p]]
    R = [[1, 0], [0, 1]]; Mt = [[1, (-Q) % p], [1, 0]]
    while k:
        if k & 1: R = mul(R, Mt)
        Mt = mul(Mt, Mt); k >>= 1
    return R[1][0]

def odd(x):
    while x % 2 == 0: x //= 2
    return x

def v2(x): return (x & -x).bit_length() - 1

def main():
    X, T0, fn = int(float(sys.argv[1])), int(sys.argv[2]), sys.argv[3]
    got = set()
    for l in open(fn):
        D, p, t, j1, j2, o1, o2 = map(int, l.split()); got.add((D, p, t, j1, j2, o1, o2))
    want = set()
    for p in primerange(3, X):
        o = n_order(2, p)
        if any(q > T0 for q in factorint(odd(o))): continue
        for D in DL:
            Q = (1 - D) // 4
            if D % p == 0 or Q % p == 0: continue
            js = jacobi_symbol(D % p, p)
            m = p - js
            z = min(d for d in divisors(m) if matU(d, Q, p) == 0)
            o1, o2 = odd(o), odd(z)
            if any(q > T0 for q in factorint(o2)): continue
            from math import gcd
            if js == 1 and gcd(o1, o2) != 1: continue
            want.add((D, p, int(js == 1), v2(o), v2(z), o1, o2))
    print(f"X = {X}, T0 = {T0}: census3 records {len(got)}, brute-force {len(want)}, "
          f"missing {len(want - got)}, extra {len(got - want)}")
    for r in sorted(want - got)[:10]: print("  missing", r)
    for r in sorted(got - want)[:10]: print("  extra", r)
    for D in DL:
        print(f"  D = {D}: {sum(1 for r in want if r[0] == D)} (split {sum(1 for r in want if r[0] == D and r[2])})")

if __name__ == "__main__":
    main()
