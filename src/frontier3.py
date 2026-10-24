"""Frontier with exponent caps (and optional split primes).
Per odd prime q <= y (q != 5): side in {A, B, X} and cap >= 1. Row p is in the pool iff every q-exponent of
odd(ord_p 2) is <= cap(q) with side A, and every q-exponent of odd(z(p)) is <= cap(q) with side B.
Cost log2|G| = 1 + sum over used q of log2(q-1) + (e_q - 1) log2 q, e_q = max exponent in the pool (as before).
Rows can be restricted to y' < y (pairs whose p -+ 1 are y'-smooth) to get the growth curve from one file.
Usage: python3 frontier3.py z pairs_file yprime iters seeds [split]"""
import sys, math, random, json
from collections import defaultdict
import numpy as np
from sympy import primerange, isprime, factorint
from pool import build, order_mod, v2, odd

def fib_mod(k, p):
    a, b = 0, 1
    for bit in bin(k)[2:]:
        c = a * ((2 * b - a) % p) % p; d = (a * a + b * b) % p
        a, b = (d, (c + d) % p) if bit == '1' else (c, d)
    return a

def zrank(p, m, fac):            # rank of apparition of Fibonacci mod p, m = p - (5/p), fac = factorint(m)
    o = m
    for q, e in fac.items():
        for _ in range(e):
            if fib_mod(o // q, p) == 0: o //= q
            else: break
    return o

def build_split(fname):
    rows = []
    for l in open(fname):
        x = int(l); p = 2 * x + 1
        if not (x % 5 == 0 or (x + 1) % 5 == 0) or p == 5 or not isprime(p): continue
        fm = factorint(p - 1); o2 = order_mod(2, p, fm); zz = zrank(p, p - 1, fm)
        rows.append(dict(p=p, cls=(v2(o2), v2(zz)), R1=factorint(odd(o2)), R2=factorint(odd(zz)), split=1))
    return rows

def prep(rows, qs):
    ix = {q: i for i, q in enumerate(qs)}; n, k = len(rows), len(qs)
    E1 = np.zeros((n, k), dtype=np.int8); E2 = np.zeros((n, k), dtype=np.int8)
    for t, r in enumerate(rows):
        for q, e in r['R1'].items(): E1[t, ix[q]] = e
        for q, e in r['R2'].items(): E2[t, ix[q]] = e
    small = np.array([r['p'] in ix for r in rows]); sp = np.array([ix.get(r['p'], 0) for r in rows])
    return E1, E2, small, sp

def evaluate(P, lq1, lq, side, cap):
    E1, E2, small, sp = P
    cA = np.where(side == 0, cap, 0); cB = np.where(side == 1, cap, 0)
    sel = (E1 <= cA).all(axis=1) & (E2 <= cB).all(axis=1) & ~(small & (side[sp] != 2))
    N = int(sel.sum())
    if N == 0: return 0, 1.0, sel
    e = E1[sel].max(axis=0).astype(np.int64) + E2[sel].max(axis=0).astype(np.int64)
    return N, 1.0 + float(((e > 0) * lq1).sum() + (np.maximum(e - 1, 0) * lq).sum()), sel

def anneal(P, lq1, lq, emax, lam, iters, seed):
    rnd = np.random.default_rng(seed); k = len(lq1)
    side = rnd.integers(0, 2, k); cap = emax.copy()
    N, lg, _ = evaluate(P, lq1, lq, side, cap); cur = N - lam * lg; best = (cur, side.copy(), cap.copy(), N, lg)
    for it in range(iters):
        T = 2.0 * (1 - it / iters) + 0.02
        i = rnd.integers(k); os, oc = side[i], cap[i]
        if rnd.random() < 0.5: side[i] = rnd.choice([s for s in (0, 1, 2) if s != os])
        else: cap[i] = min(max(1, oc + rnd.choice([-1, 1])), emax[i])
        N, lg, _ = evaluate(P, lq1, lq, side, cap); val = N - lam * lg
        if val >= cur or rnd.random() < math.exp((val - cur) / T):
            cur = val
            if val > best[0]: best = (val, side.copy(), cap.copy(), N, lg)
        else: side[i], cap[i] = os, oc
    return best

if __name__ == "__main__":
    z = int(sys.argv[1]); fn = sys.argv[2]; yp = int(sys.argv[3]); iters = int(sys.argv[4]); seeds = int(sys.argv[5])
    rows = build(z, fn) + (build_split(fn) if 'split' in sys.argv else [])
    qs = [q for q in primerange(3, yp + 1) if q != 5]
    ok = lambda r: all(q <= yp and q != 5 for q in list(r['R1']) + list(r['R2']))
    rows = [r for r in rows if ok(r)]
    bycls = defaultdict(list)
    for r in rows: bycls[r['cls']].append(r)
    c = max(bycls, key=lambda c: len(bycls[c])); R = bycls[c]; P = prep(R, qs)
    emax = np.maximum(P[0].max(axis=0), P[1].max(axis=0)).astype(np.int64); emax[emax == 0] = 1
    lq1 = np.array([math.log2(q - 1) for q in qs]); lq = np.array([math.log2(q) for q in qs])
    print(f"y'={yp} (file z={z}{', +split' if 'split' in sys.argv else ''}): class {c} size {len(R)}", flush=True)
    out = []
    for lam in (0.0, 0.05, 0.08, 0.1, 0.12, 0.15, 0.2):
        b = max((anneal(P, lq1, lq, emax, lam, iters, s) for s in range(seeds)), key=lambda t: t[0])
        print(f"  lam {lam}: log2|G| {b[4]:.1f} N {b[3]} ratio {b[3]/b[4]:.3f}", flush=True)
        out.append((lam, round(b[4], 1), b[3]))
    json.dump(dict(y=yp, cls=str(c), size=len(R), pts=out), open(f"../data/frontier3_y{yp}{'_split' if 'split' in sys.argv else ''}.json", "w"))
