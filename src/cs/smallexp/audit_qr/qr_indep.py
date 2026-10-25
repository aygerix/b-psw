#!/usr/bin/env python3
"""Independent audit of the QR-relation lemma (no code shared with zsolve.py).
Claim: for any odd prime l, (r+1)^2 vectors of F_l^r give, deterministically, a nonzero relation
sum c_i x_i = 0 with every c_i a nonzero square on its (nonempty) support.
Strict version: r+1 consecutive blocks of r+1 inputs, no push-back, so the bound (r+1)^2 is exact.
nu found deterministically: l = 3 mod 4 -> nu = -1 and one representation n0 = s0 + 1 of a nonsquare
(binary search for a square s0 with s0 + 1 nonsquare); l = 1 mod 4 -> nu taken from the data (a block whose
kernel has only square entries is already a relation)."""
import sys, random

def leg(x, l):
    x %= l
    if x == 0: return 0
    return 1 if pow(x, (l - 1) // 2, l) == 1 else -1

def issq(x, l): return leg(x, l) >= 0

def kernel(cols, l, r):
    """nonzero c with sum c_j cols_j = 0 (len(cols) > rank). Plain Gaussian elimination on the r x m matrix."""
    m = len(cols)
    A = [[cols[j][i] % l for j in range(m)] for i in range(r)]
    piv = []; row = 0
    for c in range(m):
        p = next((i for i in range(row, r) if A[i][c]), None)
        if p is None: continue
        A[row], A[p] = A[p], A[row]
        inv = pow(A[row][c], l - 2, l)
        A[row] = [(v * inv) % l for v in A[row]]
        for i in range(r):
            if i != row and A[i][c]:
                f = A[i][c]; A[i] = [(u - f * v) % l for u, v in zip(A[i], A[row])]
        piv.append((row, c)); row += 1
        if row == r: break
    pc = {c for _, c in piv}
    f = next(c for c in range(m) if c not in pc)
    t = [0] * m; t[f] = 1
    for rr, c in piv: t[c] = (-A[rr][f]) % l
    return t

def find_rep_minus1(l):
    """l = 3 mod 4: square s0 with s0+1 nonsquare, by bisection on [1, l-1] (1 square, l-1 = -1 nonsquare)."""
    lo, hi = 1, l - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if issq(mid, l): lo = mid
        else: hi = mid
    assert issq(lo, l) and not issq(lo + 1, l)
    return lo, lo + 1

def qr_relation(X, l, r):
    """X: list of (r+1)^2 vectors. Returns dict index -> coefficient."""
    B = r + 1
    assert len(X) >= B * B
    blocks = [list(range(k * B, (k + 1) * B)) for k in range(B)]
    ts = [kernel([X[i] for i in blk], l, r) for blk in blocks]
    for blk, t in zip(blocks, ts):
        assert any(t) and all(sum(t[j] * X[blk[j]][i] for j in range(B)) % l == 0 for i in range(r))
    if l % 4 == 3:
        nu = l - 1; s0, n0 = find_rep_minus1(l)
    else:
        nu = None
        for blk, t in zip(blocks, ts):
            if all(issq(v, l) for v in t):
                return {blk[j]: t[j] for j in range(B) if t[j]}, 'block'
            if nu is None:
                nu = next(v for v in t if not issq(v, l))
    assert not issq(nu, l)
    dec = []
    for blk, t in zip(blocks, ts):
        a, b = {}, {}
        for j, tj in enumerate(t):
            if tj == 0: continue
            if issq(tj, l): aj, bj = tj, 0
            elif l % 4 == 1: aj, bj = 0, (-tj * pow(nu, l - 2, l)) % l
            else:                                   # t = a + b (nu = -1), t nonsquare: t = (t/n0) s0 + (t/n0)
                q = tj * pow(n0, l - 2, l) % l; aj, bj = q * s0 % l, q
            assert issq(aj, l) and issq(bj, l) and (aj - nu * bj - tj) % l == 0
            if aj: a[blk[j]] = aj
            if bj: b[blk[j]] = bj
        w = [sum(c * X[i][k] for i, c in b.items()) % l for k in range(r)]
        aw = [sum(c * X[i][k] for i, c in a.items()) % l for k in range(r)]
        assert all((aw[k] - nu * w[k]) % l == 0 for k in range(r))
        dec.append((w, a, b))
    for w, a, b in dec:
        if not any(w):
            return (b if b else a), 'zeroblock'
    lam = kernel([w for w, _, _ in dec], l, r)
    rel = {}
    for (w, a, b), lk in zip(dec, lam):
        if lk == 0: continue
        if issq(lk, l):
            for i, c in b.items(): rel[i] = lk * c % l
        else:
            s = lk * pow(nu, l - 2, l) % l
            for i, c in a.items(): rel[i] = s * c % l
    return rel, 'LA'

def check(rel, X, l, r):
    assert rel, "empty relation"
    assert all(c % l and issq(c, l) for c in rel.values()), "non-square coefficient"
    tot = [sum(c * X[i][k] for i, c in rel.items()) % l for k in range(r)]
    assert not any(tot), "relation does not vanish"

def gens(l, r, rng):
    N = (r + 1) ** 2
    if l < 10 ** 4:
        nsq = [x for x in range(1, l) if not issq(x, l)]
    else:
        nsq = []
        while len(nsq) < 50:
            x = rng.randrange(1, l)
            if not issq(x, l): nsq.append(x)
    yield 'uniform', [[rng.randrange(l) for _ in range(r)] for _ in range(N)]
    d = max(1, r // 3)
    basis = [[rng.randrange(l) for _ in range(r)] for _ in range(d)]
    yield 'rank-r/3', [[sum(rng.randrange(l) * bb[k] for bb in basis) % l for k in range(r)] for _ in range(N)]
    yield 'zeros', [[0] * r for _ in range(N)]
    v = [rng.randrange(1, l) for _ in range(r)]
    yield 'repeated', [list(v) for _ in range(N)]
    yield 'sparse', [[(rng.randrange(1, l) if rng.random() < 2.0 / r else 0) for _ in range(r)] for _ in range(N)]
    X = []                                           # blocks: basis + nonsquare combination -> nonsquare kernels
    for k in range(r + 1):
        P = [[rng.randrange(l) for _ in range(r)] for _ in range(r)]
        # random invertible-ish basis; last vector = sum nonsquare * basis
        coeffs = [rng.choice(nsq) for _ in range(r)]
        last = [sum(coeffs[j] * P[j][i] for j in range(r)) % l for i in range(r)]
        X += P + [last]
    yield 'nonsq-kernels', X
    X = []                                           # all blocks give the same w direction; scaled by nonsquares
    e = [1] + [0] * (r - 1)
    for k in range(r + 1):
        s = rng.choice(nsq)
        X += [[(s * v) % l for v in e]] + [[0] * r for _ in range(r)]
    yield 'collinear', X
    yield 'unit-vectors', [[1 if k == (i % r) else 0 for k in range(r)] for i in range(N)]

if __name__ == '__main__':
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
    primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 37, 41, 43, 1000003, 998244353, 1000000007]
    stats = {}
    for l in primes:
        for r in [1, 2, 3, 5, 8, 13, 20]:
            if l > 10 ** 5 and r > 13: continue
            for _ in range(3):
                for name, X in gens(l, r, rng):
                    rel, how = qr_relation(X, l, r)
                    check(rel, X, l, r)
                    stats[(l, how)] = stats.get((l, how), 0) + 1
    for k in sorted(stats): print(k, stats[k])
    print("ALL CHECKS PASSED")
