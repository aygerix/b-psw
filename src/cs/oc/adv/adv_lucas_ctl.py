#!/usr/bin/env python3
"""Positive control for the Lucas-half criterion used in adv_union.py (per-prime j2 and odd z from adv_rows.lucas_side):
all products of two distinct primes < PMAX; criterion (one j2, odd z | n + 1, j2 <= v2(n+1), with D = Selfridge D of n)
compared with src/bpsw.py strong_lucas. Usage: adv_lucas_ctl.py PMAX"""
import sys, os, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bpsw import strong_lucas, selfridge_D
from adv_rows import mul, pw, ONE, jac
from sympy import primerange, factorint
def v2(x): return (x & -x).bit_length() - 1
cache = {}
def zdata(p, D):
    key = (p, D)
    if key in cache: return cache[key]
    Q = (1 - D) // 4
    if Q % p == 0 or D % p == 0: cache[key] = None; return None
    eps = jac(D, p); N = p - eps
    i2 = (p + 1) // 2
    g = mul((i2, i2), (i2, i2), D, p); iq = pow(Q % p, -1, p); g = (g[0] * iq % p, g[1] * iq % p)
    z = N
    for r, e in factorint(N).items():
        for _ in range(e):
            if pw(g, z // r, D, p) == ONE: z //= r
            else: break
    cache[key] = (v2(z), z >> v2(z)); return cache[key]
P = list(primerange(3, int(sys.argv[1])))
tot = agree = pos = 0
for a, b in itertools.combinations(P, 2):
    n = a * b; D = selfridge_D(n)
    if D is None: continue
    za, zb = zdata(a, D), zdata(b, D)
    if za is None or zb is None: crit = False
    else: crit = za[0] == zb[0] and za[0] <= v2(n + 1) and (n + 1) % za[1] == 0 and (n + 1) % zb[1] == 0
    sl = strong_lucas(n); tot += 1; agree += (sl == crit); pos += sl
    if sl != crit: print('DISAGREE', a, b, D, crit, sl)
print(f'pairs {tot}, criterion agrees {agree}, strong Lucas psp found {pos}')
