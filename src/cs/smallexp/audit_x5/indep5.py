#!/usr/bin/env python3
"""Independent audit of the exponent-5 claim (parts A-C do not import zsolve).

A: the norm-trick step in F_5: every t in F_5 is a - 2b with a, b in {0, +-1}, a*b = 0 (exhaustive).
B: independent implementation of  ternary relations -> doubling -> doubling -> LA  for F_5^r WITHOUT
   recycling: consumes exactly (r+1)^5 base elements (unless a degenerate early exit gives a zero-sum
   sooner).  Checks: indices distinct, exact sum 0 mod 5, every copy of every k-fold element has the
   claimed value, ternary relations have coefficients in {+-1}.
C: adversarial inputs (many repeated vectors / low-rank pools / pools containing zero vectors) for B.
Usage: python3 indep5.py rmax trials seed
"""
import sys, random, time
import numpy as np

P = 5

def kernel(cols):
    """cols: m vectors in F_5^r (list of np arrays).  Return nonzero c in F_5^m with sum c_j col_j = 0 (r < m)."""
    M = np.array(cols, dtype=np.int64).T % P            # r x m
    r, m = M.shape
    M = M.copy(); pivcols = []; row = 0
    for c in range(m):
        if row >= r: break
        nz = [i for i in range(row, r) if M[i, c] % P]
        if not nz: continue
        i = nz[0]
        M[[row, i]] = M[[i, row]]
        inv = pow(int(M[row, c]), P - 2, P)
        M[row] = (M[row] * inv) % P
        for i2 in range(r):
            if i2 != row and M[i2, c]:
                M[i2] = (M[i2] - M[i2, c] * M[row]) % P
        pivcols.append(c); row += 1
    free = [c for c in range(m) if c not in pivcols]
    assert free, "no kernel"
    f = free[-1]                                        # deliberately a different choice than zsolve (last free)
    x = np.zeros(m, dtype=np.int64); x[f] = 1
    for k, pc in enumerate(pivcols):
        x[pc] = (-M[k, f]) % P
    assert not ((np.array(cols, dtype=np.int64).T @ x) % P).any()
    return x

SPLIT = {1: (1, 0), 4: (4, 0), 2: (0, 4), 3: (0, 1)}   # t = a - 2b (mod 5)

def partA():
    for t in range(5):
        if t == 0: continue
        a, b = SPLIT[t]
        assert (a - 2 * b) % P == t and a in (0, 1, 4) and b in (0, 1, 4) and a * b == 0
    print("A: F_5 split t = a - 2b, a,b in {0,+-1}, disjoint: OK")

class Pool:
    def __init__(self, X):
        self.X, self.i = X, 0
    def take(self):
        i = self.i; self.i += 1
        return i

class ZeroSumFound(Exception):
    def __init__(self, idx): self.idx = idx

class KElt:
    """k-fold element: copies = list of disjoint index tuples, each summing to value v (mod 5)."""
    def __init__(self, v, copies): self.v, self.copies = v, copies

def csum(X, idx):
    if not len(idx): return np.zeros(X.shape[1], dtype=np.int64)
    return X[np.array(idx, dtype=np.int64)].sum(axis=0) % P

def ternary(src, r, X):
    """(r+1) blocks of (r+1) elements from src(); returns (Kplus, Kminus) lists of elements with
    sum_{K+} v = sum_{K-} v, K+ u K- nonempty, all coefficients +-1."""
    blocks = []
    for _ in range(r + 1):
        E = [src() for _ in range(r + 1)]
        t = kernel([e.v for e in E])
        A = [(e, SPLIT[int(ti)][0]) for e, ti in zip(E, t) if ti and SPLIT[int(ti)][0]]
        B = [(e, SPLIT[int(ti)][1]) for e, ti in zip(E, t) if ti and SPLIT[int(ti)][1]]
        w = sum((c * e.v for e, c in B), np.zeros(r, dtype=np.int64)) % P
        wa = sum((c * e.v for e, c in A), np.zeros(r, dtype=np.int64)) % P
        assert not ((wa - 2 * w) % P).any()             # sum a x = 2 * sum b x
        if not w.any():                                 # zero block: a relation by itself
            rel = B if B else A
            return [e for e, c in rel if c == 1], [e for e, c in rel if c == 4]
        blocks.append((w, A, B))
    lam = kernel([w for w, _, _ in blocks])
    rel = []
    for (w, A, B), lk in zip(blocks, lam):
        lk = int(lk)
        if lk == 0: continue
        if lk in (1, 4):
            rel += [(e, (lk * c) % P) for e, c in B]
        else:
            s = (lk * 3) % P                            # lk / 2 = lk * 3 mod 5
            rel += [(e, (s * c) % P) for e, c in A]
    assert rel and all(c in (1, 4) for _, c in rel)
    tot = sum((c * e.v for e, c in rel), np.zeros(r, dtype=np.int64)) % P
    assert not tot.any()
    return [e for e, c in rel if c == 1], [e for e, c in rel if c == 4]

def double(src, r, X, m):
    """One 2m-fold element from one ternary relation over m-fold elements."""
    Kp, Km = ternary(src, r, X)
    if not Km or not Kp:                                # one-sided relation: 0/1 zero-sum of m-fold elements
        K = Kp or Km
        idx = tuple(i for e in K for i in e.copies[0])
        raise ZeroSumFound(idx)
    W = sum((e.v for e in Kp), np.zeros(r, dtype=np.int64)) % P
    copies = []
    for K in (Kp, Km):
        for j in range(m):
            copies.append(tuple(i for e in K for i in e.copies[j]))
    for c in copies:
        assert not ((csum(X, c) - W) % P).any(), "copy value mismatch"
    allidx = [i for c in copies for i in c]
    assert len(allidx) == len(set(allidx)), "copies not disjoint"
    return KElt(W, copies)

def solve(X, r):
    pool = Pool(X)
    def base():
        i = pool.take()
        return KElt(X[i] % P, [(i,)])
    try:
        two = lambda: double(base, r, X, 1)
        four = lambda: double(two, r, X, 2)
        F = [four() for _ in range(r + 1)]
        lam = kernel([e.v for e in F])
        idx = []
        for e, lk in zip(F, lam):
            for j in range(int(lk)):
                idx.extend(e.copies[j])
        idx = tuple(idx)
    except ZeroSumFound as z:
        idx = z.idx
    assert idx and len(idx) == len(set(idx)), "not a 0/1 subset"
    assert not csum(X, idx).any(), "sum nonzero"
    return pool.i, len(idx)

def main():
    rmax, trials, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    rng = random.Random(seed); g = np.random.default_rng(seed)
    partA()
    for r in range(1, rmax + 1):
        B = (r + 1) ** 5
        for kind in ("random", "lowrank", "repeats", "zeros"):
            t0 = time.time(); used = []; sizes = []
            for _ in range(trials):
                if kind == "random":
                    X = g.integers(0, 5, size=(B, r), dtype=np.int64)
                elif kind == "lowrank":                 # pool spans a subspace of rank r-1
                    Bas = g.integers(0, 5, size=(max(1, r - 1), r), dtype=np.int64)
                    C = g.integers(0, 5, size=(B, Bas.shape[0]), dtype=np.int64)
                    X = (C @ Bas) % P
                elif kind == "repeats":                 # only 3 distinct nonzero vectors, no two summing simply
                    V = [np.array([rng.randrange(1, 5) for _ in range(r)]) for _ in range(3)]
                    X = np.array([V[rng.randrange(3)] for _ in range(B)], dtype=np.int64)
                else:
                    X = g.integers(0, 5, size=(B, r), dtype=np.int64)
                    X[g.random(B) < 0.05] = 0
                u, s = solve(X, r)
                assert u <= B
                used.append(u); sizes.append(s)
            print(f"B r={r:2d} {kind:8s} bound={B:8d} consumed max={max(used):8d} mean={sum(used)/len(used):10.1f}"
                  f" |I| mean={sum(sizes)/len(sizes):7.1f} [{time.time()-t0:.1f}s]", flush=True)

if __name__ == "__main__":
    main()
