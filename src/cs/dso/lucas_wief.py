"""Lucas-Wieferich check (independent of src/wief/lw.gp): for the base-2 Wieferich primes p = 1093, 3511 and every
Selfridge D (5, -7, 9, -11, ...) with |D| <= DMAX, is z_D(p^2) = z_D(p)?  (U_n(1, Q), Q = (1 - D)/4.)
p^2 | n in a BPSW psp needs this for the selected D (and p Wieferich).
Usage: python3 lucas_wief.py DMAX
"""
import sys
from math import isqrt
DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1000

def U(k, Q, m):
    """U_k(1, Q) mod m by 2x2 matrix power [[1, -Q], [1, 0]]^k, entry (2,1)."""
    a, b, c, d = 1, 0, 0, 1          # identity
    x, y, z, w = 1, (-Q) % m, 1, 0   # M
    while k:
        if k & 1:
            a, b, c, d = (a * x + b * z) % m, (a * y + b * w) % m, (c * x + d * z) % m, (c * y + d * w) % m
        x, y, z, w = (x * x + y * z) % m, (x * y + y * w) % m, (z * x + w * z) % m, (z * y + w * w) % m
        k >>= 1
    return c

def jac(a, n):
    a %= n; r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: r = -r
        a %= n
    return r if n == 1 else 0

def factor(n):
    f = {}; q = 2
    while q * q <= n:
        while n % q == 0: f[q] = f.get(q, 0) + 1; n //= q
        q += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def entry(p, D, Q):
    m = p - jac(D, p)
    for q in factor(m):
        while m % q == 0 and U(m // q, Q, p) == 0:
            m //= q
    assert U(m, Q, p) == 0
    return m

def _useq(Q, n):
    u = [0, 1]
    for _ in range(n): u.append(u[-1] - Q * u[-2])
    return u
for _Q in (-1, 2, 3, -3, 4):
    _u = _useq(_Q, 40)
    assert all(U(k, _Q, 10**30) == _u[k] % 10**30 for k in range(40)), 'U self-test failed'
for p in (1093, 3511):
    assert pow(2, p - 1, p * p) == 1, "not Wieferich"
    hits = []; tested = 0; skipped = []
    D = 5
    while abs(D) <= DMAX:
        Q = (1 - D) // 4
        sq = D > 0 and isqrt(D) ** 2 == D
        if D % p == 0 or Q % p == 0:
            skipped.append(D)
        else:
            tested += 1
            z = entry(p, D, Q)
            if U(z, Q, p * p) == 0:
                hits.append((D, z, jac(D, p), 'square D (never selected)' if sq else ''))
        D = -D - 2 if D > 0 else -D + 2
    print(f"p = {p}: ord_p(2) = {next(k for k in range(1, p) if pow(2, k, p) == 1)}, Selfridge D with |D| <= {DMAX}: tested {tested}, "
          f"skipped (p | D or p | Q) {skipped}; Lucas-Wieferich hits [D, z(p), (D/p)]: {hits}")
