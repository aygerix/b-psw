#!/usr/bin/env python3
"""src/cs/e3/adv/v5_scan.py -- independent re-scan of the E3 rho-family with src/bpsw.py (no gp code reused).

Family (rho = a/b):  B = 2x^2 - (4 rho + 1) x + 2 rho - 2,  f1 = 1 + B,  f2 = 1 + B x/(2 rho),  f3 = 1 + B (x-1)/(2 rho - 1).
x = X/d.  The X-range replicates rho_scan.gp's loop exactly: M = lcm of the three denominators,
X in [M*floor(-Xmax/M), M*floor(-Xmax/M) + M*K + M - 1], K = floor((Xmax - lo)/M).
Keeps triples of integer values >= 5, pairwise distinct, all prime (sympy.isprime, deterministic < 2^64).
For each: n = p1 p2 p3, bpsw.selfridge_D, bpsw.strong_fermat2, bpsw.strong_lucas, Korselt, and on all-inert triples
the Lucas defects dL_j = odd(z_j)/gcd(odd(z_j), n+1), dLt = lcm_j odd(z_j) / gcd(., n+1)  (z = rank of apparition).

Usage: v5_scan.py SRC AMAX BMAX DLIST XMAX OUT [--estimate]
"""
import sys, math, time
from fractions import Fraction as Fr
import numpy as np

SRC, AMAX, BMAX, DLIST, XMAX, OUT = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), [int(t) for t in sys.argv[4].split(',')], int(float(sys.argv[5])), sys.argv[6]
EST = '--estimate' in sys.argv
sys.path.insert(0, SRC)
import bpsw
from sympy import isprime, factorint

def pmul(p, q):
    r = [Fr(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            r[i + j] += x * y
    return r
def padd(p, q):
    n = max(len(p), len(q)); return [(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)]
def pscale(p, s): return [s * c for c in p]

def family(rho):
    B = [2 * rho - 2, -(4 * rho + 1), Fr(2)]
    f1 = padd([Fr(1)], B)
    f2 = padd([Fr(1)], pscale(pmul(B, [Fr(0), Fr(1)]), 1 / (2 * rho)))
    f3 = padd([Fr(1)], pscale(pmul(B, [Fr(-1), Fr(1)]), 1 / (2 * rho - 1)))
    return [f1, f2, f3]

def intpoly(f, d):
    g = [c / Fr(d) ** k for k, c in enumerate(f)]
    L = 1
    for c in g: L = L * c.denominator // math.gcd(L, c.denominator)
    return [int(c * L) for c in g], L

def okmask(num, L):
    r = np.arange(L, dtype=np.int64); acc = np.zeros(L, dtype=np.int64)
    for c in reversed(num):
        acc = (acc * r + (c % L)) % L
    return acc == 0

def ev(num, L, X):
    v = 0
    for c in reversed(num): v = v * X + c
    q, rem = divmod(v, L); assert rem == 0
    return q

def oddp(m):
    while m % 2 == 0: m //= 2
    return m

def zrank(p, Q, D):
    N = p - bpsw.jacobi(D % p, p); z = N
    for l, e in factorint(N).items():
        for _ in range(e):
            if bpsw.lucas_UV(z // l, 1, Q, p)[0] == 0: z //= l
            else: break
    return z

t0 = time.time()
fams = 0; integral = 0; rows = []
for b in range(1, BMAX + 1):
    for a in range(-AMAX, AMAX + 1):
        if a == 0 or math.gcd(a, b) != 1 or 2 * a == b: continue
        f = family(Fr(a, b))
        for d in DLIST:
            fams += 1
            P = [intpoly(fj, d) for fj in f]
            M = 1
            for _, L in P: M = M * L // math.gcd(M, L)
            lo = M * ((-XMAX) // M); K = (XMAX - lo) // M; hi = lo + M * K + M - 1
            X = np.arange(lo, hi + 1, dtype=np.int64)
            mask = np.ones(len(X), dtype=bool)
            for num, L in P: mask &= okmask(num, L)[X % L]
            Xs = X[mask]; integral += len(Xs)
            if EST: continue
            for Xv in Xs.tolist():
                v = [ev(num, L, Xv) for num, L in P]
                if min(v) < 5 or len(set(v)) < 3: continue
                if not all(isprime(t) for t in v): continue
                rows.append((a, b, d, Xv, v))
print(f"families x d: {fams}  integral X: {integral}  prime triples: {len(rows)}  time {time.time() - t0:.1f}s", flush=True)
if EST: sys.exit(0)
stats = dict(tri=0, inert=0, s2=0, sl=0, bp=0, kor=0)
lr = []; lr64 = []; dl1max = 0; dl1ratio = 0.0; dl1one = 0; exph = 0.0; best = None; dl23 = []
with open(OUT, 'w') as fh:
    for a, b, d, Xv, v in rows:
        n = v[0] * v[1] * v[2]
        assert n == v[0] * v[1] * v[2] and all(isprime(t) for t in v)
        D = bpsw.selfridge_D(n); Q = (1 - D) // 4
        s2 = bpsw.strong_fermat2(n); sl = bpsw.strong_lucas(n)
        kor = all((n - 1) % (t - 1) == 0 for t in v)
        stats['tri'] += 1; stats['s2'] += s2; stats['sl'] += sl; stats['bp'] += (s2 and sl); stats['kor'] += kor
        if s2 and sl: print("BPSW PASS", n, v, a, b, d, Xv, flush=True)
        inert = all(bpsw.jacobi(D % t, t) == -1 for t in v)
        rec = [n] + v + [f"{a}/{b}", d, Xv, D, int(s2), int(sl), int(kor)]
        if inert:
            stats['inert'] += 1
            z = [oddp(zrank(t, Q, D)) for t in v]
            dl = [zj // math.gcd(zj, n + 1) for zj in z]
            Lz = 1
            for zj in z: Lz = Lz * zj // math.gcd(Lz, zj)
            dLt = Lz // math.gcd(Lz, n + 1)
            r = math.log(dLt) / math.log(n) if dLt > 1 else 0.0
            lr.append(r)
            if n > 2 ** 64: lr64.append(r)
            dl1max = max(dl1max, dl[0]); dl1one += (dl[0] == 1)
            if dl[0] > 1: dl1ratio = max(dl1ratio, math.log(dl[0]) / math.log(v[0]))
            dl23.append((math.log(max(dl[1], 1)) / math.log(v[1]), math.log(max(dl[2], 1)) / math.log(v[2])))
            exph += 1.0 / dLt
            if best is None or dLt < best[0]: best = (dLt, n, v, a, b, d, Xv, D)
            rec += dl + [dLt]
        fh.write(' '.join(map(str, rec)) + '\n')
print(stats)
lr.sort(); lr64.sort()
print(f"inert {len(lr)}: log dLt/log n min {lr[0]:.3f} median {lr[len(lr)//2]:.3f} max {lr[-1]:.3f}")
if lr64: print(f"inert with n > 2^64 ({len(lr64)}): min {lr64[0]:.3f} median {lr64[len(lr64)//2]:.3f}")
print(f"member-1 dL1: max {dl1max}, max log dL1/log p1 {dl1ratio:.3f}, dL1 == 1 on {dl1one}/{len(lr)} inert")
print(f"members 2,3 on inert: mean log dL2/log p2 {sum(t[0] for t in dl23)/len(dl23):.4f}  mean log dL3/log p3 {sum(t[1] for t in dl23)/len(dl23):.4f}")
print(f"sum over inert of 1/dLt = {exph:.4e};  smallest dLt: {best}")
print(f"time {time.time() - t0:.1f}s")
