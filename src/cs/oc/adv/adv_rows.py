#!/usr/bin/env python3
"""Adversarial recomputation of oc_fam.gp rows with new code.
Row format: fam par bits D eps j1 j2 oo200 oz200 oo1500 oz1500 roughF roughL.
Members regenerated from (fam, par) by independent formulas; Lucas side in F_p[r]/(r^2 - D), gamma = alpha^2/Q.
Usage: adv_rows.py MAXBITS out_file in_file [in_file ...]   (rows with bits > MAXBITS are skipped)"""
import sys, math, time
from sympy import divisors, primerange, mobius as _mob

def mob(n): return int(_mob(n))
def v2(x): return (x & -x).bit_length() - 1
def fibpair(n):  # (F_n, F_{n+1})
    if n == 0: return (0, 1)
    a, b = fibpair(n >> 1)
    c = a * (2 * b - a); d = a * a + b * b
    return (d, c + d) if n & 1 else (c, d)
def F(n): return fibpair(n)[0]
def L(n): a, b = fibpair(n); return 2 * b - a
def member(fam, par):
    if fam == 'mer': return 2**par - 1
    if fam == 'wag': return (2**par + 1) // 3
    if fam == 'fib': return F(par)
    if fam == 'luc': return L(par)
    if fam in ('gm+', 'gm-'):
        k = (par - 1) // 2; return 2**par + (1 if fam == 'gm+' else -1) * 2**(k + 1) + 1
    if fam == 'cyc':
        num = den = 1
        for d in divisors(par):
            mu = mob(par // d)
            if mu == 1: num *= 2**d - 1
            elif mu == -1: den *= 2**d - 1
        P = num // den; assert P * den == num; return P // math.gcd(P, par)
    if fam == 'fcy':
        num = den = 1
        for d in divisors(par):
            mu = mob(par // d)
            if mu == 1: num *= F(d)
            elif mu == -1: den *= F(d)
        P = num // den; assert P * den == num; return P // math.gcd(P, par)
    if fam.startswith('rie'): return int(fam[3:]) * 2**par - 1
    if fam.startswith('pro'): return int(fam[3:]) * 2**par + 1
    raise ValueError(fam)

PR = {B: list(primerange(3, B + 1)) for B in (200, 1500)}
def split_smooth(m, B):
    S = 1
    for r in PR[B]:
        while m % r == 0: m //= r; S *= r
    return S, m

def mul(a, b, D, p): return ((a[0]*b[0] + D*a[1]*b[1]) % p, (a[0]*b[1] + a[1]*b[0]) % p)
def pw(a, e, D, p):
    r = (1, 0)
    while e:
        if e & 1: r = mul(r, a, D, p)
        a = mul(a, a, D, p); e >>= 1
    return r
ONE = (1, 0)

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

def fermat_side(p):
    N = p - 1; v = v2(N); m = N >> v
    y = pow(2, m, p); j1 = 0
    while y != 1: y = y * y % p; j1 += 1
    h = pow(2, 1 << v, p)          # odd-order part
    out = {}
    for B in (200, 1500):
        S, R = split_smooth(m, B)
        if pow(h, S, p) != 1: out[B] = 0
        else:
            o = S
            for r in PR[B]:
                if o % r: continue
                while o % r == 0 and pow(h, o // r, p) == 1: o //= r
            out[B] = o
    S, R = split_smooth(m, 1500)
    rough = 0.0 if R == 1 or pow(2, N // R, p) == 1 else math.log2(R)
    return j1, out[200], out[1500], rough

def lucas_side(p, D):
    Q = (1 - D) // 4
    eps = jac(D, p); N = p - eps; v = v2(N); m = N >> v
    i2 = (p + 1) // 2
    g = mul((i2, i2), (i2, i2), D, p); iq = pow(Q % p, -1, p); g = (g[0] * iq % p, g[1] * iq % p)
    y = pw(g, m, D, p); j2 = 0
    while y != ONE:
        y = mul(y, y, D, p); j2 += 1
        if j2 > v: raise RuntimeError
    h = g
    for _ in range(v): h = mul(h, h, D, p)
    out = {}
    for B in (200, 1500):
        S, R = split_smooth(m, B)
        if pw(h, S, D, p) != ONE: out[B] = 0
        else:
            o = S
            for r in PR[B]:
                if o % r: continue
                while o % r == 0 and pw(h, o // r, D, p) == ONE: o //= r
            out[B] = o
    S, R = split_smooth(m, 1500)
    rough = 0.0 if R == 1 or pw(g, N // R, D, p) == ONE else math.log2(R)
    return eps, j2, out[200], out[1500], rough

if __name__ == '__main__':
    MAXB = int(sys.argv[1]); outf = open(sys.argv[2], 'w'); t0 = time.time()
    nrow = nok = nbad = nskip = 0
    cache = {}
    for path in sys.argv[3:]:
        for line in open(path):
            t = line.split()
            if len(t) != 13: continue
            fam, par, bits, D = t[0], int(t[1]), int(t[2]), int(t[3])
            if bits > MAXB: nskip += 1; continue
            key = (fam, par)
            if key not in cache:
                p = member(fam, par)
                cache.clear(); cache[key] = (p, fermat_side(p))
            p, fs = cache[key]
            if p.bit_length() != bits: print('BITS MISMATCH', fam, par, bits, p.bit_length()); nbad += 1; continue
            Q = (1 - D) // 4
            if Q % p == 0:
                print('DEGENERATE p | Q', fam, par, D, 'their row', t[4:]); continue
            eps, j2, oz200, oz1500, rL = lucas_side(p, D)
            j1, oo200, oo1500, rF = fs
            mine = [eps, j1, j2, oo200, oz200, oo1500, oz1500]
            theirs = list(map(int, t[4:11]))
            okr = (mine == theirs) and abs(rF - float(t[11])) < 0.06 and abs(rL - float(t[12])) < 0.06
            nrow += 1
            if okr: nok += 1
            else: nbad += 1; print('ROW MISMATCH', fam, par, D, 'theirs', t[4:], 'mine', mine, round(rF, 1), round(rL, 1))
            outf.write(f'{fam} {par} {bits} {D} {eps} {j1} {j2} {oo200} {oz200} {oo1500} {oz1500} {rF:.1f} {rL:.1f}\n')
    print(f'rows checked {nrow}, agree {nok}, mismatch {nbad}, skipped (> {MAXB} bits) {nskip}, {time.time()-t0:.1f}s')
