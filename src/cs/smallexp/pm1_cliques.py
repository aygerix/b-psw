#!/usr/bin/env python3
"""+-1-closed subpools: subsets C of one (D, class) block such that for every member p and every row prime power
q^e required by any member (q^e | odd(ord) [F side] or odd(z) [L side]), p == +-1 mod q^e, and sides are
consistent.  For such C the generated group has exponent 2, so the BPSW congruences for any S in C become
GF(2)-linear (sufficient form over the rows of C):  F-row: #{p in S: p == -1} even;  L-row: #{p in S: p == +1} even;
|S| odd.  We find maximal cliques of the pairwise compatibility graph (Bron-Kerbosch), report |C|, the number of
rows, the GF(2) rank, and solve the system; any solution n is checked with src/bpsw.py.
Usage: python3 pm1_cliques.py export_file [D ...]
"""
import sys, os, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
import numpy as np
import bpsw
import zsolve as z
from pool_group import parse, qprime

def rows_of(B, i):
    Q, r = B['Q'], B['P'][i]
    return [(qprime(Q[qi]), e, 'F') for qi, e in r['F']] + [(qprime(Q[qi]), e, 'L') for qi, e in r['L']]

def pm1(p, q, e):
    m = q ** e; x = p % m
    return x == 1 or x == m - 1

def nullspace2(A):
    A = A.copy() % 2; rows, cols = A.shape
    piv = z._rref(A, 2); pc = [c for _, c in piv]
    basis = []
    for f in range(cols):
        if f in pc: continue
        v = np.zeros(cols, dtype=np.int64); v[f] = 1
        for r, c in piv: v[c] = A[r, f] % 2
        basis.append(v)
    return basis

def cliques(adj, n):
    """Bron-Kerbosch with pivoting over bitmasks; yields maximal cliques."""
    out = []
    def bk(R, P, X):
        if not P and not X:
            out.append(R); return
        u = (P | X).bit_length() - 1
        cand = P & ~adj[u]
        while cand:
            v = cand.bit_length() - 1; cand &= ~(1 << v)
            bk(R | 1 << v, P & adj[v], X & adj[v])
            P &= ~(1 << v); X |= 1 << v
    bk(0, (1 << n) - 1, 0)
    return out

def main():
    Ds = set(map(int, sys.argv[2:])) if len(sys.argv) > 2 else None
    best = []
    for B in parse(sys.argv[1]):
        if Ds and B['D'] not in Ds: continue
        P = B['P']; n = len(P)
        if n < 3: continue
        R = [rows_of(B, i) for i in range(n)]
        adj = [0] * n
        for i in range(n):
            for j in range(n):
                if i == j or j in P[i]['conf']: continue
                if all(pm1(P[j]['p'], q, e) for q, e, s in R[i]) and all(pm1(P[i]['p'], q, e) for q, e, s in R[j]) \
                   and P[i]['p'] not in {q for q, _, _ in R[j]} and P[j]['p'] not in {q for q, _, _ in R[i]}:
                    adj[i] |= 1 << j
        cl = cliques(adj, n)
        for c in cl:
            S = [i for i in range(n) if c >> i & 1]
            if len(S) < 3: continue
            rows = sorted({(q, e, s) for i in S for (q, e, s) in R[i]})
            # merge prime powers: keep the max exponent per (q, side)
            mx = {}
            for q, e, s in rows: mx[(q, s)] = max(mx.get((q, s), 0), e)
            rws = sorted(mx.items())
            M = []
            for (q, s), e in rws:
                M.append([1 if (P[i]['p'] % q ** e == (q ** e - 1 if s == 'F' else 1)) else 0 for i in S])
            M.append([1] * len(S))
            tgt = [0] * len(rws) + [1]
            A = np.array(M, dtype=np.int64) % 2
            rank = len(z._rref(A.copy(), 2))
            sol = z.affine_solve([np.array(col) for col in np.array(M).T], np.array(tgt), 2)
            found = None
            if sol is not None:
                ker = nullspace2(A)
                for mask in range(1 << min(len(ker), 16)):
                    x = sol.copy() % 2
                    for b in range(len(ker)):
                        if mask >> b & 1: x = (x + ker[b]) % 2
                    T = [S[j] for j in range(len(S)) if x[j]]
                    if len(T) >= 3:
                        nn = math.prod(P[i]['p'] for i in T); found = (len(T), bpsw.bpsw_pass(nn), nn)
                        if found[1]: break
            best.append((len(S) - rank, len(S), len(rws), rank, B['D'], B['cls'], found))
            if found:
                print(f"SOLUTION D={B['D']} class {B['cls']} |C|={len(S)} rows={len(rws)} : |S|={found[0]} bpsw={found[1]} n={found[2]}", flush=True)
    best.sort(key=lambda t: (-t[0], -t[1]))
    print("top cliques by slack |C| - rank (slack >= 1 with the parity row means a GF(2) solution exists):")
    for t in best[:15]:
        print(f"  slack={t[0]:3d} |C|={t[1]:3d} rows={t[2]:3d} rank={t[3]:3d} D={t[4]} class={t[5]} sol={t[6] and t[6][:2]}")
    sizes = [t[1] for t in best]
    print(f"maximal cliques with |C|>=3: {len(best)}, max |C| = {max(sizes) if sizes else 0}, "
          f"max slack = {best[0][0] if best else None}")

if __name__ == "__main__":
    main()
