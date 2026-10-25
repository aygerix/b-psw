#!/usr/bin/env python3
"""Group structure of designed-row BPSW pools and an end-to-end solver run.

Input: an export file of src/c3/export.py (one block per (D, 2-adic class)).  For each class:
  - maximum independent set of the side-conflict graph (exact branch and bound),
  - the group G_S = (Z/A_S)^* x (Z/B_S)^* x Z/2 (parity) of that set: bits, exponent, layer ranks, the layered
    bound prod f_l(rho) (finite only if exp(G_S) has primes in {2,3,5}),
  - heuristic expected number of hits  E = sum over independent sets S (odd #inert, |S| >= 2) of 1/|G_S|,
  - solver: zsolve target (forced parity element) on the max independent set when exp(G) is {2,3,5}-smooth,
    else exhaustive subset check; every candidate n is run through the verbatim checker src/bpsw.py.
Usage: python3 pool_group.py export_file [D ...]
"""
import sys, os, math, itertools
from collections import defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from sympy import primitive_root, factorint, totient
from sympy.ntheory import discrete_log
import bpsw
import zsolve as z

def parse(path):
    blocks, f = [], open(path)
    for line in f:
        t = line.split()
        if not t: continue
        if t[0] == 'CLASS':
            D, j1, j2, k, nq = map(int, t[1:6]); Q = list(map(int, f.readline().split()[1:]))
            P = []
            for _ in range(k):
                s = list(map(int, f.readline().split()[1:]))
                p, eps = s[0], s[1]; pos = 2 + nq
                nF = s[pos]; Fr = [(s[pos + 1 + 2 * i], s[pos + 2 + 2 * i]) for i in range(nF)]; pos += 1 + 2 * nF
                nL = s[pos]; Lr = [(s[pos + 1 + 2 * i], s[pos + 2 + 2 * i]) for i in range(nL)]; pos += 1 + 2 * nL
                nc = s[pos]; conf = set(s[pos + 1: pos + 1 + nc])
                P.append(dict(p=p, eps=eps, F=Fr, L=Lr, conf=conf))
            blocks.append(dict(D=D, cls=(j1, j2), Q=Q, P=P))
    return blocks

def qprime(m):
    f = factorint(m); assert len(f) == 1; return next(iter(f))

def max_indep(P):
    n = len(P); adj = [0] * n
    for i, r in enumerate(P):
        for j in r['conf']: adj[i] |= 1 << j
    best = [0, 0]
    def rec(cand, cur, size):
        if size + bin(cand).count('1') <= best[0]: return
        if not cand:
            best[0], best[1] = size, cur; return
        i = cand.bit_length() - 1
        rec(cand & ~adj[i] & ~(1 << i), cur | 1 << i, size + 1)
        rec(cand & ~(1 << i), cur, size)
    rec((1 << n) - 1, 0, 0)
    return [i for i in range(n) if best[1] >> i & 1]

def indep_sets(P):
    n = len(P); adj = [0] * n
    for i, r in enumerate(P):
        for j in r['conf']: adj[i] |= 1 << j
    def rec(start, cur, forb):
        yield cur
        for i in range(start, n):
            if not (forb >> i & 1):
                yield from rec(i + 1, cur + [i], forb | adj[i])
    yield from rec(0, [], 0)

def group_of(B, S):
    """Row moduli (q^e, side) for a subset S of the class block B."""
    Q = B['Q']; need = {}
    for i in S:
        for side, req in (('F', B['P'][i]['F']), ('L', B['P'][i]['L'])):
            for qi, e in req:
                q = qprime(Q[qi]); key = (q, side)
                need[key] = max(need.get(key, 0), e)
    return need

def log2_order(need):
    return sum(math.log2(int(totient(q ** e))) for (q, s), e in need.items()) + 1.0

def analyse(B, run_solver=True):
    P = B['P']; k = len(P)
    S = max_indep(P)
    need = group_of(B, S)
    rowq = {q for q, _ in need}
    S = [i for i in S if P[i]['p'] not in rowq]          # a pool prime that is itself a row cannot divide n
    need = group_of(B, S)
    mods = [int(totient(q ** e)) for (q, s), e in sorted(need.items())] + [2]
    G = z.Group(mods)
    expo = math.lcm(*mods)
    ok235 = set(factorint(expo)) <= {2, 3, 5}
    bnd = z.bound(G) if ok235 else float('inf')
    lay = [(l, kk, len(a)) for l, kk, a in G.layers()]
    E, nsets = 0.0, 0
    if k <= 40:
        for T in indep_sets(P):
            ninert = sum(1 for i in T if P[i]['eps'] == -1)
            if len(T) >= 2 and ninert % 2 == 1:
                E += 2.0 ** (-log2_order(group_of(B, T))); nsets += 1
    print(f"D={B['D']} class {B['cls']}: k={k} max-indep={len(S)} rows={len(need)} log2|G_S|={math.log2(math.prod(mods)):.1f} "
          f"exp={expo} ({'2,3,5-smooth' if ok235 else 'other primes'}) layers={lay[:8]}{'...' if len(lay) > 8 else ''} "
          f"bound={bnd if bnd < 1e30 else 'inf'}  odd indep sets={nsets} E[hits]={E:.3g}", flush=True)
    hits = []
    if run_solver and len(S) >= 3:
        # elements x_p: dlogs of p (F rows) and -p (L rows), plus parity 1; forced u = parity 1
        keys = sorted(need)
        gens = {key: primitive_root(key[0] ** need[key]) for key in keys}
        X = []
        for i in S:
            p = P[i]['p']; v = []
            for key in keys:
                m = key[0] ** need[key]; a = (p if key[1] == 'F' else -p) % m
                v.append(discrete_log(m, a, gens[key]))
            X.append(tuple(v) + (1,))
        u = tuple([0] * len(keys) + [1])
        cand = []
        if ok235:
            J = z.Solver(G, X).target(u)
            print(f"   solver: {'no solution from the pool (exhausted)' if J is None else 'subset ' + str(J)}")
            if J: cand.append([S[j] for j in J])
        if len(S) <= 22:                                  # exhaustive cross-check inside the max independent set
            for r in range(3, len(S) + 1, 2):
                for T in itertools.combinations(range(len(S)), r):
                    tot = [0] * len(mods)
                    for j in T: tot = [(a + b) % m for a, b, m in zip(tot, X[j], mods)]
                    if all(t == 0 for t in tot[:-1]) and tot[-1] == 1:
                        cand.append([S[j] for j in T])
        for T in cand:
            n = math.prod(P[i]['p'] for i in T)
            ok = bpsw.bpsw_pass(n)
            print(f"   candidate n = {n} ({len(T)} primes): strong_fermat2 = {bpsw.strong_fermat2(n)}, bpsw_pass = {ok}")
            if ok: hits.append(n)
    return hits

def main():
    Ds = set(map(int, sys.argv[2:])) if len(sys.argv) > 2 else None
    allhits = []
    for B in parse(sys.argv[1]):
        if Ds and B['D'] not in Ds: continue
        if len(B['P']) < 2: continue
        allhits += analyse(B)
    print("HITS:", allhits)

if __name__ == "__main__":
    main()
