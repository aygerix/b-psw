"""Vectorised version of frontier.py (same objective, same annealer) for large pair files.
Rows -> bitmasks over odd primes q <= z (q != 5) and exponent matrices; pool membership is
(mask1 & ~A) == 0 and (mask2 & ~B) == 0.  Usage: python3 frontier_np.py z pairs_file [iters]"""
import sys, math, random, json
from collections import defaultdict
import numpy as np
from sympy import primerange
from pool import build

def prep(rows, qs):
    ix = {q: i for i, q in enumerate(qs)}
    n = len(rows); k = len(qs)
    m1 = np.zeros(n, dtype=np.uint64); m2 = np.zeros(n, dtype=np.uint64)
    E1 = np.zeros((n, k), dtype=np.int8); E2 = np.zeros((n, k), dtype=np.int8)
    ps = np.array([r['p'] for r in rows], dtype=object)
    for t, r in enumerate(rows):
        for q, e in r['R1'].items():
            m1[t] |= np.uint64(1 << ix[q]); E1[t, ix[q]] = e
        for q, e in r['R2'].items():
            m2[t] |= np.uint64(1 << ix[q]); E2[t, ix[q]] = e
    small = [t for t, r in enumerate(rows) if r['p'] in ix]          # pool primes that are also modulus primes
    return m1, m2, E1, E2, small

def evaluate(P, qs, lq1, lq, assign):
    m1, m2, E1, E2, small = P
    A = sum(1 << i for i, q in enumerate(qs) if assign[q] == 'A')
    B = sum(1 << i for i, q in enumerate(qs) if assign[q] == 'B')
    sel = ((m1 & np.uint64(~A & (2**64 - 1))) == 0) & ((m2 & np.uint64(~B & (2**64 - 1))) == 0)
    for t in small:
        if assign[qs_p[t]] != 'X': sel[t] = False
    N = int(sel.sum())
    if N == 0: return 0, 1.0
    e1 = E1[sel].max(axis=0); e2 = E2[sel].max(axis=0)
    e = e1.astype(np.int64) + e2.astype(np.int64)        # a prime sits on one side only
    lg = 1.0 + float(((e > 0) * lq1).sum() + (np.maximum(e - 1, 0) * lq).sum())
    return N, lg

def anneal_l(P, qs, lq1, lq, lam, iters, seed):
    rnd = random.Random(seed)
    assign = {q: rnd.choice('AB') for q in qs}
    N, lg = evaluate(P, qs, lq1, lq, assign); cur = N - lam * lg
    best = (cur, dict(assign), N, lg)
    for it in range(iters):
        T = 2.0 * (1 - it / iters) + 0.02
        q = rnd.choice(qs); old = assign[q]
        assign[q] = rnd.choice([c for c in 'ABX' if c != old])
        N, lg = evaluate(P, qs, lq1, lq, assign); val = N - lam * lg
        if val >= cur or rnd.random() < math.exp((val - cur) / T):
            cur = val
            if val > best[0]: best = (val, dict(assign), N, lg)
        else: assign[q] = old
    return best

if __name__ == "__main__":
    z = int(sys.argv[1]); rows = build(z, sys.argv[2]); iters = int(sys.argv[3]) if len(sys.argv) > 3 else 4000
    qs = [q for q in primerange(3, z + 1) if q != 5]
    lq1 = np.array([math.log2(q - 1) for q in qs]); lq = np.array([math.log2(q) for q in qs])
    bycls = defaultdict(list)
    for r in rows: bycls[r['cls']].append(r)
    print(f"z={z}: inert rows {len(rows)}", {str(c): len(v) for c, v in sorted(bycls.items(), key=lambda t: -len(t[1]))[:6]}, flush=True)
    res = {}
    for c in sorted(bycls, key=lambda c: -len(bycls[c]))[:3]:
        R = bycls[c]; P = prep(R, qs); qs_p = {t: R[t]['p'] for t in P[4]}
        pts = []
        for lam in (0.0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75):
            b = max((anneal_l(P, qs, lq1, lq, lam, iters, s) for s in range(3)), key=lambda t: t[0])
            pts.append((round(b[3], 1), b[2], round(b[2] - b[3], 1)))
            print(f"  class {c} lam {lam}: log2|G| {b[3]:.1f} N {b[2]}", flush=True)
        print(f"z={z} class {c} (|class|={len(R)}) frontier (log2|G|, N, N-log2|G|):", sorted(set(pts)), flush=True)
        res[str(c)] = sorted(set(pts))
    json.dump(res, open(f"../data/frontier_np_z{z}.json", "w"))
