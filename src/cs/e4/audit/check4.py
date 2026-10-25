#!/usr/bin/env python3
"""Independent check of the four 'usable' forced primes p > 2^33 named in the E4 result, the Mersenne classes,
the mod-4-rule exclusion of 486077598803, and two cheap products through src/bpsw.py."""
import sys
sys.path.insert(0, sys.argv[1] if len(sys.argv) > 1 else '../../..')
from bpsw import bpsw_pass
from sympy import factorint, isprime, jacobi_symbol


def fibmod(n, m):
    def fp(n):
        if n == 0:
            return (0, 1)
        a, b = fp(n >> 1)
        c = a * ((2 * b - a) % m) % m
        d = (a * a + b * b) % m
        return (d, (c + d) % m) if n & 1 else (c, d)
    return fp(n)[0]


def order(N, ok):
    n = N
    for q, a in factorint(N).items():
        for _ in range(a):
            if ok(n // q):
                n //= q
            else:
                break
    return n


def v2(n):
    return (n & -n).bit_length() - 1


def fs(n):
    return '*'.join('%d^%d' % (q, a) if a > 1 else str(q) for q, a in sorted(factorint(n).items()))


cases = [('M107', 2 ** 107 - 1), ('M127', 2 ** 127 - 1), ('LF167', 18104700793), ('LF1596', 486077598803)]
for name, p in cases:
    assert isprime(p)
    o = order(p - 1, lambda n: pow(2, n, p) == 1)
    z = order(p + 1, lambda n: fibmod(n, p) == 0)
    print('%s p=%d bits=%d (5/p)=%d p%%8=%d' % (name, p, p.bit_length(), jacobi_symbol(5, p), p % 8))
    print('   ord2=%s  z=%s  class=(%d,%d)  max odd prime: ord2 %s, z %s' % (
        fs(o), fs(z), v2(o), v2(z), max([q for q in factorint(o) if q > 2] or [1]), max([q for q in factorint(z) if q > 2] or [1])))
    print('   p-1=%s   p+1=%s' % (fs(p - 1), fs(p + 1)))
    print('   bpsw_pass(p) =', bpsw_pass(p))
# mod-4 (Pomerance) side rule for class (1,2): p == 3 mod 8, A=(p-1)/2 all primes == 1 mod 4 (F side),
# odd part of the Lucas order (B = (p+1)/4) all primes == 3 mod 4. With true orders: odd ord2 primes == 1 mod 4,
# odd z primes == 3 mod 4.
p = 486077598803
o = order(p - 1, lambda n: pow(2, n, p) == 1)
z = order(p + 1, lambda n: fibmod(n, p) == 0)
print('mod-4 rule check for', p)
print('   odd ord2 primes and residues mod 4:', [(q, q % 4) for q in factorint(o) if q > 2])
print('   odd z primes and residues mod 4:   ', [(q, q % 4) for q in factorint(z) if q > 2])
# sanity: products of different-class primes are not BPSW psp
M107, M127 = 2 ** 107 - 1, 2 ** 127 - 1
print('bpsw_pass(M107*M127) =', bpsw_pass(M107 * M127))
print('bpsw_pass(18104700793*486077598803) =', bpsw_pass(18104700793 * 486077598803))
# is z(M_q) = 2^q exactly? (F_{2^(q-1)} != 0 mod M_q)
for q in (107, 127):
    Mq = 2 ** q - 1
    print('q=%d: F_{2^q} == 0 mod M_q: %s ; F_{2^(q-1)} == 0 mod M_q: %s' % (q, fibmod(2 ** q, Mq) == 0, fibmod(2 ** (q - 1), Mq) == 0))
