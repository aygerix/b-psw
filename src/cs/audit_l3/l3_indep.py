#!/usr/bin/env python3
"""Independent re-implementation (no code shared with smallexp/zsolve.py) of the exponent-3 doubled-element
zero-sum solver in (Z/3)^r.  Two modes:
  norecycle : blocks are disjoint consecutive runs of r+1 elements; at most (r+1)^2 elements are read.
  recycle   : elements with kernel coefficient 0 go back into the pool.
Every output is verified: nonempty, distinct indices, exact sum 0 in (Z/3)^r.
Inputs: random and adversarial (structured) sequences.  Usage: l3_indep.py <seed>"""
import sys, random

def kernel_vec(vecs, r):
    """vecs: list of m tuples in F_3^r. Return nonzero c in F_3^m with sum c_j vecs_j = 0, or None.
    Gaussian elimination tracking combinations (bitless, plain lists)."""
    m = len(vecs)
    rows = [list(v) + [1 if k == j else 0 for k in range(m)] for j, v in enumerate(vecs)]
    piv_rows = []          # (pivot col, row)
    for row in rows:
        row = row[:]
        for pc, pr in piv_rows:
            if row[pc]:
                f = row[pc] * pr[pc] % 3     # pr[pc] = 1 after normalisation, keep general
                row = [(x - f * y) % 3 for x, y in zip(row, pr)]
        lead = next((i for i in range(r) if row[i]), None)
        if lead is None:
            return row[r:]
        inv = row[lead]  # 1 or 2, self-inverse mod 3
        row = [(x * inv) % 3 for x in row]
        piv_rows.append((lead, row))
    return None

def add(x, y):
    return tuple((a + b) % 3 for a, b in zip(x, y))

def vsum(X, I, r):
    s = (0,) * r
    for i in I:
        s = add(s, X[i])
    return s

def solve(X, r, recycle):
    """Return (I, elements read) or (None, read)."""
    pool = []           # recycled indices
    nxt = [0]
    def take(k):
        out = []
        while len(out) < k:
            if pool:
                out.append(pool.pop())
            elif nxt[0] < len(X):
                out.append(nxt[0]); nxt[0] += 1
            else:
                return None
        return out
    blocks = []         # (w, S2, S12)
    while len(blocks) < r + 1:
        E = take(r + 1)
        if E is None:
            return None, nxt[0]
        c = kernel_vec([X[i] for i in E], r)
        assert c is not None and any(c)
        S1 = [i for i, ci in zip(E, c) if ci == 1]
        S2 = [i for i, ci in zip(E, c) if ci == 2]
        if recycle:
            pool.extend(i for i, ci in zip(E, c) if ci == 0)
        # identity: sum(S1) + 2 sum(S2) = 0  =>  sum(S1) = sum(S2)
        assert vsum(X, S1, r) == vsum(X, S2, r)
        w = vsum(X, S2, r)
        if not any(w):
            I = S2 if S2 else S1
            return sorted(I), nxt[0]
        blocks.append((w, S2, S1 + S2))
    lam = kernel_vec([b[0] for b in blocks], r)
    assert lam is not None and any(lam)
    I = []
    for (w, S2, S12), lj in zip(blocks, lam):
        if lj == 1:
            I += S2
        elif lj == 2:
            I += S12
    return sorted(I), nxt[0]

def check(X, r, I):
    assert I and len(set(I)) == len(I), "empty or repeated"
    assert vsum(X, I, r) == (0,) * r, "nonzero sum"

def gen(kind, r, n, rng):
    e = lambda i, a=1: tuple(a if k == i else 0 for k in range(r))
    if kind == 'random':
        return [tuple(rng.randrange(3) for _ in range(r)) for _ in range(n)]
    if kind == 'const':
        v = tuple(rng.randrange(1, 3) if k == 0 else rng.randrange(3) for k in range(r))
        return [v] * n
    if kind == 'basis_cycle':                       # e1,e1,e2,e2,...,er,er repeated (longest zs-free prefix 2r)
        return [e((j // 2) % r) for j in range(n)]
    if kind == 'basis_rand':
        return [e(rng.randrange(r), rng.randrange(1, 3)) for _ in range(n)]
    if kind == 'lowrank':                           # inside a random rank-2 subspace
        a = tuple(rng.randrange(3) for _ in range(r)); b = tuple(rng.randrange(3) for _ in range(r))
        out = []
        for _ in range(n):
            s, t = rng.randrange(3), rng.randrange(3)
            out.append(tuple((s * x + t * y) % 3 for x, y in zip(a, b)))
        return out
    if kind == 'sparse_last':                       # mostly zero rows; nonzero only in last coordinate block
        return [tuple(rng.randrange(3) if k >= r - 2 else 0 for k in range(r)) for _ in range(n)]
    if kind == 'with_zeros':
        return [((0,) * r) if rng.random() < 0.3 else tuple(rng.randrange(3) for _ in range(r)) for _ in range(n)]
    if kind == 'twos':                              # all entries in {0,2}
        return [tuple(2 * rng.randrange(2) for _ in range(r)) for _ in range(n)]
    raise ValueError(kind)

KINDS = ['random', 'const', 'basis_cycle', 'basis_rand', 'lowrank', 'sparse_last', 'with_zeros', 'twos']

def main():
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
    # 1. 'always': norecycle mode on exactly (r+1)^2 elements, every input family
    for r in (1, 2, 3, 5, 8, 13, 20):
        n = (r + 1) ** 2
        for kind in KINDS:
            for _ in range(200 if r <= 8 else 40):
                X = gen(kind, r, n, rng)
                for rec in (False, True):
                    I, used = solve(X, r, rec)
                    assert I is not None, (r, kind, rec)
                    assert used <= n
                    check(X, r, I)
        print(f"r={r:3d}: all {len(KINDS)} families x both modes solved within (r+1)^2={n}", flush=True)
    # 2. the recycled count on random inputs (reproduce 273 / 2487 / 6845)
    for r, tr in ((20, 40), (60, 8), (100, 3)):
        us, us0 = [], []
        for _ in range(tr):
            X = gen('random', r, (r + 1) ** 2, rng)
            I, used = solve(X, r, True); check(X, r, I); us.append(used)
            I, used = solve(X, r, False); check(X, r, I); us0.append(used)
        m = sum(us) / len(us)
        print(f"r={r:3d} random: recycle mean {m:8.1f} max {max(us)}  ratio to (r+1)^2 {m/(r+1)**2:.3f};"
              f" norecycle mean {sum(us0)/len(us0):8.1f}", flush=True)

if __name__ == '__main__':
    main()
