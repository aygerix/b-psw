#!/usr/bin/env python3
"""Independent audit of the layered GF(2) zero-sum algorithm for abelian 2-groups G = prod Z/2^{a_j}.
Written from scratch (does not import zsolve).  Bitmask GF(2) elimination.
Modes:
  norecycle : each layer consumes fixed blocks of r_i+1 lower items; unused items discarded.
              Uses exactly prod_i (r_i+1) base elements; must ALWAYS succeed (deterministic bound).
  recycle   : unused items of a kernel go back to the lower stream (as zsolve does).
Inputs: random, and an adversarial family for (Z/4)^r designed so every kernel has full support.
Every output is checked: distinct indices, exact sum 0 in G."""
import sys, random, time

def dep(vecs):
    """vecs: list of int bitmasks. Return a nonempty subset (list of positions) with XOR 0, using the
    first linearly dependent position (smallest prefix), or None."""
    basis = {}  # pivot bit -> (vec, combo bitmask over positions)
    for p, v in enumerate(vecs):
        c = 1 << p
        while v:
            h = v.bit_length() - 1
            if h in basis:
                bv, bc = basis[h]; v ^= bv; c ^= bc
            else:
                basis[h] = (v, c); break
        if v == 0:
            return [q for q in range(p + 1) if (c >> q) & 1]
    return None

class Solver:
    def __init__(self, a, X, recycle):
        self.a, self.X, self.recycle = a, X, recycle
        self.K = max(a)
        self.ranks = [sum(1 for e in a if e >= i) for i in range(1, self.K + 1)]
        self.nxt = 0
        self.buf = [[] for _ in range(self.K + 1)]
    def img(self, val, i):          # layer-i image (i = 1..K) of val in 2^{i-1} G
        m = 0; t = 0
        for j, e in enumerate(self.a):
            if e >= i:
                assert val[j] % (1 << (i - 1)) == 0
                if (val[j] >> (i - 1)) & 1: m |= 1 << t
                t += 1
            else:
                assert val[j] % (1 << e) == 0   # automatically 0
        return m
    def pull(self, lev):
        """An item (idxs, val) with val in 2^lev G."""
        if self.buf[lev]:
            return self.buf[lev].pop()
        if lev == 0:
            if self.nxt >= len(self.X): return None
            i = self.nxt; self.nxt += 1
            return ((i,), self.X[i])
        i = lev; r = self.ranks[i - 1]
        items = []
        for _ in range(r + 1):
            it = self.pull(lev - 1)
            if it is None: return None
            items.append(it)
            # (no early exit: block of exactly r+1, as in the deterministic bound)
        S = dep([self.img(v, i) for _, v in items])
        assert S is not None and len(S) > 0       # r+1 vectors in F_2^r are dependent
        idx = []; val = [0] * len(self.a)
        for q in S:
            idx.extend(items[q][0])
            for j, e in enumerate(self.a): val[j] = (val[j] + items[q][1][j]) % (1 << e)
        if self.recycle:
            Sset = set(S)
            for q in reversed(range(len(items))):
                if q not in Sset: self.buf[lev - 1].append(items[q])
        val = tuple(val)
        for j, e in enumerate(self.a): assert val[j] % (1 << min(i, e)) == 0   # val in 2^i G
        return (tuple(idx), val)
    def zero_sum(self):
        return self.pull(self.K)

def check(a, X, I):
    assert I and len(set(I)) == len(I)
    for j, e in enumerate(a):
        assert sum(X[i][j] for i in I) % (1 << e) == 0, "sum not zero"

def adversarial_Z4(r, rng):
    """(Z/4)^r: blocks of r+1 elements with mod-2 images e_1..e_r, (1..1); lifts chosen so the t-th block's
    layer-2 image is e_t (t < r+1) or (1..1) (t = r+1).  Every kernel has full support for the 'first
    dependent prefix' rule, so both modes consume (r+1)^2."""
    X = []
    for t in range(r + 1):
        blk = []
        for s in range(r + 1):
            v = [0] * r
            if s < r: v[s] = 1
            else: v = [1] * r
            blk.append(v)
        # integer sum of mod-2 parts per coordinate = 2, so carry = (1..1); add 2*b on first element
        target = [0] * r
        if t < r: target[t] = 1
        else: target = [1] * r
        b = [(1 + target[j]) % 2 for j in range(r)]
        blk[0] = [(blk[0][j] + 2 * b[j]) % 4 for j in range(r)]
        X.extend(tuple(v) for v in blk)
    return X

def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    rng = random.Random(seed)
    cases = [([2] * 30, 50), ([4] * 5, 200), ([4] * 20, 50), ([4] * 60, 10), ([4] * 100, 3),
             ([8] * 4, 200), ([8] * 10, 20), ([16] * 3, 100), ([2, 4, 8], 300), ([2, 2, 4, 8, 8, 16], 50),
             ([32] * 2, 100)]
    for a, trials in cases:
        ranks = [sum(1 for e in [x.bit_length() - 1 for x in a] if e >= i) for i in range(1, max(x.bit_length() - 1 for x in a) + 1)]
        ex = [x.bit_length() - 1 for x in a]
        B = 1
        for r in ranks: B *= r + 1
        t0 = time.time()
        res = {}
        for mode in ("norecycle", "recycle"):
            used = []
            for _ in range(trials):
                X = [tuple(rng.randrange(1 << e) for e in ex) for _ in range(B)]
                S = Solver(ex, X, mode == "recycle")
                out = S.zero_sum()
                assert out is not None, f"FAILED with {B} elements"   # bound must suffice
                check(ex, X, out[0])
                used.append(S.nxt)
                if mode == "norecycle": assert S.nxt == B
            res[mode] = (sum(used) / len(used), max(used))
        r0 = ranks[0]
        print(f"a={sorted(set(a))}x{len(a)} ranks={ranks} bound={B} norecycle used={res['norecycle'][0]:.0f} "
              f"recycle mean={res['recycle'][0]:.1f} max={res['recycle'][1]} ratio_to_bound={res['recycle'][0]/B:.3f}  {time.time()-t0:.1f}s",
              flush=True)
    for r in (5, 20, 60):
        X = adversarial_Z4(r, rng)
        for mode in (False, True):
            S = Solver([2] * r, X, mode)
            out = S.zero_sum()
            assert out is not None
            check([2] * r, X, out[0])
            print(f"adversarial (Z/4)^{r} recycle={mode}: used {S.nxt} of (r+1)^2={(r+1)**2}, |I|={len(out[0])}", flush=True)

if __name__ == "__main__":
    main()
