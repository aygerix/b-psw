#!/usr/bin/env python3
"""Adversarial check of Lemma OC4 and Lemma OC3 with new code (F_p[sqrt D] arithmetic, no Lucas U/V formulas).
OC4: p inert, p == 3 mod 4, (Q/p) = 1, s^2 = Q: v2(z) = v2(p+1) - 1  <=>  ((2Q+s)/p) = -1; symbol sign-free.
OC3: p inert => min(v2 ord_p 2, v2 z(p)) <= 1.
Also checks Lemma 18(b) on the (Q/p) = -1 inert primes (v2(z) = v2(p+1)).
Usage: adv_oc4.py PMAX NRAND BITS"""
import sys, random, time
from sympy import primerange, isprime

def mul(a, b, D, p): return ((a[0]*b[0] + D*a[1]*b[1]) % p, (a[0]*b[1] + a[1]*b[0]) % p)
def pw(a, e, D, p):
    r = (1, 0)
    while e:
        if e & 1: r = mul(r, a, D, p)
        a = mul(a, a, D, p); e >>= 1
    return r
def v2(x): return (x & -x).bit_length() - 1
def leg(a, p):
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1

def v2z_inert(p, D, Q):
    i2 = (p + 1) // 2
    al = (i2, i2); g = mul(al, al, D, p); iq = pow(Q % p, -1, p); g = (g[0]*iq % p, g[1]*iq % p)
    N = p + 1; g = pw(g, N >> v2(N), D, p); j = 0
    while g != (1, 0):
        g = mul(g, g, D, p); j += 1
        if j > v2(N): raise RuntimeError
    return j
def v2ord2(p):
    m = (p - 1) >> v2(p - 1); y = pow(2, m, p); j = 0
    while y != 1: y = y * y % p; j += 1
    return j

DL = [5,-7,-11,13,-15,17,-19,-23,29,-31,37,41,-43,-47,53,-55,-59,61,-67,-71,73]
stats = {D: [0, 0, 0, 0, 0, 0] for D in DL}  # inert, oc3 viol, oc4 tested, oc4 bad, sign-dep, l18 bad

def test(p, D):
    Q = (1 - D) // 4
    if (2 * D * Q) % p == 0: return
    if leg(D, p) != -1: return
    st = stats[D]; st[0] += 1
    j2 = v2z_inert(p, D, Q); j1 = v2ord2(p)
    if min(j1, j2) > 1: st[1] += 1; print('OC3 VIOLATION', D, p, j1, j2)
    k = v2(p + 1)
    if leg(Q, p) == -1:
        if j2 != k: st[5] += 1; print('L18 VIOLATION', D, p)
        return
    if p % 4 != 3: return
    s = pow(Q % p, (p + 1) // 4, p); assert s * s % p == Q % p
    a = leg(2*Q + s, p); b = leg(2*Q - s, p)
    st[2] += 1
    if a != b: st[4] += 1
    if (j2 == k - 1) != (a == -1): st[3] += 1; print('OC4 MISMATCH', D, p, j2, k, a)

if __name__ == '__main__':
    PMAX, NR, BITS = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    t0 = time.time()
    for p in primerange(3, PMAX):
        for D in DL: test(p, D)
    print(f'small primes < {PMAX} done {time.time()-t0:.1f}s', flush=True)
    random.seed(20261025)
    for D in DL:
        c = 0
        while c < NR:
            p = random.getrandbits(BITS) | (1 << (BITS - 1)) | 3   # p == 3 mod 4 bias to hit OC4
            if not isprime(p): continue
            c += 1; test(p, D)
    print(f'random {BITS}-bit primes done {time.time()-t0:.1f}s')
    tot = [0]*6
    for D in DL:
        st = stats[D]; tot = [x + y for x, y in zip(tot, st)]
        print(f'D={D:4d} inert {st[0]:7d} OC3viol {st[1]} | OC4 tested {st[2]:6d} bad {st[3]} signdep {st[4]} | L18 bad {st[5]}')
    print('TOTAL inert', tot[0], 'OC3viol', tot[1], 'OC4 tested', tot[2], 'OC4 bad', tot[3], 'signdep', tot[4], 'L18 bad', tot[5])
