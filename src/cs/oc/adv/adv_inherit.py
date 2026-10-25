#!/usr/bin/env python3
"""Independent inheritance control (line OC claim: Fibonacci primitive parts inherit +14..+62 smooth bits on the free side,
'under 10% of log2 p'). For fcy/cyc members with 2^20 < p < 2^BMAX: log2 of the 1500-smooth part of odd(p -+ 1),
against KC random primes of the same bit length in the same residue class (cyc: 1 mod e; fcy: eps mod f, eps = (5/p)).
Usage: adv_inherit.py BMAX KC"""
import sys, os, math, random, statistics
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from adv_rows import member, jac
from sympy import isprime, primerange
PR = list(primerange(3, 1501))
def smbits(N):
    while N % 2 == 0: N //= 2
    s = 0.0
    for r in PR:
        while N % r == 0: N //= r; s += math.log2(r)
    return s
BMAX, KC = int(sys.argv[1]), int(sys.argv[2])
random.seed(777)
def ctl(b, modu, res):
    cs = cm = 0.0
    for _ in range(KC):
        while True:
            r = random.getrandbits(b) | (1 << (b - 1)); r = r - (r % modu) + res
            if r < (1 << (b - 1)): r += modu
            if r.bit_length() == b and isprime(r): break
        cs += smbits(r + 1); cm += smbits(r - 1)
    return cs / KC, cm / KC
out = []
for fam, rng in (('cyc', range(4, 2001)), ('fcy', range(4, 3001))):
    for e in rng:
        if isprime(e): continue
        p = member(fam, e)
        if p < 2**20 or p.bit_length() > BMAX or not isprime(p): continue
        eps = 1 if fam == 'cyc' else jac(5, p)
        res = 1 if fam == 'cyc' else eps % e
        c = ctl(p.bit_length(), e, res)
        out.append((fam, e, p.bit_length(), eps, smbits(p + 1) - c[0], smbits(p - 1) - c[1]))
for fam in ('cyc', 'fcy'):
    for lo, hi in ((0, 64), (64, 128), (128, 256), (256, 512), (512, 1024)):
        B = [r for r in out if r[0] == fam and lo <= r[2] < hi]
        if not B: continue
        mb = statistics.mean(r[2] for r in B)
        print(f'{fam} [{lo},{hi}) n={len(B):3d} split={sum(1 for r in B if r[3] == 1):3d}  excess p+1 {statistics.mean(r[4] for r in B):+6.1f}  excess p-1 {statistics.mean(r[5] for r in B):+6.1f}  bits {mb:.0f}  p-1 excess/bits {statistics.mean(r[5] for r in B)/mb:.3f}')
