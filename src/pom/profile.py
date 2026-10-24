"""Pomerance pool (D = 5, class (1,2)) from pomgen: true orders and greedy frontier.

For each pool prime p (p = 4B - 1, A = 2B - 1, A and B T-smooth with the mod-4 rule) compute
o = ord_p(2) and z = z(p) (Fibonacci rank), check class (v2 o, v2 z) = (1, 2) and the side rule
(odd o composed of primes == 1 mod 4, odd z of primes == 3 mod 4).  The subset-product group of a
sub-pool S is G(S) = prod_q (Z/q^e_q)^*, e_q = max exponent of q in odd(o) (F side) or odd(z) (L side) over S;
the 2-adic part (n == 3 mod 8) is automatic for |S| odd.  log2 E(#solutions) ~ |S| - 1 - log2|G(S)|.
Greedy: repeatedly drop the (prime, exponent-level) whose removal most increases |S| - log2|G(S)|.
Usage: python3 profile.py poolfile
"""
import sys, math
from collections import defaultdict
from sympy import factorint

def v2(x): return (x & -x).bit_length() - 1

def order_mod(a, p, fac):
    o = p - 1
    for q, e in fac.items():
        for _ in range(e):
            if pow(a, o // q, p) == 1: o //= q
            else: break
    return o

def fib_mod(k, p):
    # F_k mod p by fast doubling
    def fd(n):
        if n == 0: return (0, 1)
        a, b = fd(n >> 1)
        c = a * ((2 * b - a) % p) % p
        d = (a * a + b * b) % p
        return (d, (c + d) % p) if n & 1 else (c, d)
    return fd(k)[0]

def zrank(p, fac):
    z = p + 1
    for q, e in fac.items():
        for _ in range(e):
            if fib_mod(z // q, p) == 0: z //= q
            else: break
    return z

def load(fn):
    ps = sorted(int(l) for l in open(fn) if l.strip())
    rows = []
    for p in ps:
        if p < 5: continue
        fm, fp = factorint(p - 1), factorint(p + 1)
        o, z = order_mod(2, p, fm), zrank(p, fp)
        if (v2(o), v2(z)) != (1, 2): continue
        F = factorint(o >> v2(o)); L = factorint(z >> v2(z))
        assert all(q % 4 == 1 for q in F) and all(q % 4 == 3 for q in L)
        rows.append((p, F, L))
    return rows

def bits(q, e): return math.log2(q - 1) + (e - 1) * math.log2(q)

def logG(rows):
    em = {}
    for _, F, L in rows:
        for q, e in list(F.items()) + list(L.items()):
            em[q] = max(em.get(q, 0), e)
    return sum(bits(q, e) for q, e in em.items()), em

def main():
    rows = load(sys.argv[1])
    B, em = logG(rows)
    idx = sorted(rows, key=lambda r: r[0])
    print(f"pool {len(rows)} primes in class (1,2), true-order log2|G| = {B:.1f}, moduli {len(em)}, "
          f"N - 1 - log2|G| = {len(rows) - 1 - B:.1f}")
    cur = list(rows); best = (len(cur) - 1 - B, len(cur), B); curve = [(B, len(cur))]
    while True:
        B, em = logG(cur)
        # candidate removals: for each q, drop all primes attaining the max exponent of q
        cand = []
        for q, e in em.items():
            drop = [r for r in cur if r[1].get(q, 0) == e or r[2].get(q, 0) == e]
            rest = [r for r in cur if not (r[1].get(q, 0) == e or r[2].get(q, 0) == e)]
            B2, _ = logG(rest)
            cand.append(((len(rest) - B2) - (len(cur) - B), q, e, len(drop), rest, B2))
        cand.sort(key=lambda c: -c[0])
        if not cand or not cur: break
        g, q, e, nd, rest, B2 = cand[0]
        cur = rest; curve.append((B2, len(cur)))
        if len(cur) - 1 - B2 > best[0]: best = (len(cur) - 1 - B2, len(cur), B2)
    pts = [c for c in curve if c[0] >= 40]
    mr = max(pts, key=lambda c: c[1] / c[0]) if pts else (0, 0)
    print("frontier (log2|G|, N):", " ".join(f"({b:.0f},{n})" for b, n in curve[::max(1, len(curve)//12)]))
    print(f"max N/log2|G| on frontier (bits >= 40): {mr[1]/max(mr[0],1):.3f} at ({mr[0]:.0f} bits, N = {mr[1]})")
    B, em = logG(cur)
    print(f"greedy frontier end: N = {len(cur)}, log2|G| = {B:.1f}, N - 1 - log2|G| = {len(cur) - 1 - B:.1f}; "
          f"best seen {best[0]:.1f} at (N, bits) = ({best[1]}, {best[2]:.1f})")

if __name__ == "__main__":
    main()
