#!/usr/bin/env python3
"""Polynomial-time 0/1 zero-sum (and target) solver for finite abelian groups whose exponent has
only the primes 2, 3, 5 (any depth), by layering over the l-adic filtrations.

Group: G = Z/m_1 + ... + Z/m_s.  Base elements: vectors x_i in G.  Output: a nonempty set I of indices
with sum_{i in I} x_i = 0 (or, with a forced element u, sum_I x_i + u = 0).

Per layer (prime l, depth k, rank rho) the layer quotient is F_l^rho; one OUTPUT is a subset of the
input elements whose sum vanishes in that quotient:
  l = 2 : GF(2) kernel of rho+1 inputs.                                   f_2(rho) = rho + 1
  l = 3 : 'norm trick' (QR-relations): blocks of rho+1 inputs, kernel t, write t_i = a_i - nu b_i with
          a_i, b_i squares; block value w = sum b_i x_i, and sum a_i x_i = nu w, so the block realises
          every multiple of w; LA over rho+1 blocks.  For l = 3 the squares are {0,1}.  f_3 = (rho+1)^2
  l = 5 : the same QR solver gives TERNARY relations (squares of F_5 are {0,+-1}) from (rho+1)^2 inputs;
          a ternary relation K+ / K- doubles multiplicity (m-fold -> 2m-fold); two doublings give
          4-fold elements, which realise every coefficient 0..4; LA over rho+1 of them.  f_5 = (rho+1)^5
A layer's outputs are the next layer's inputs, so the deterministic count is prod_layers f_l(rho).
Unused inputs are recycled (pushed back), so the realised count is smaller.
"""
import numpy as np

# ---------------------------------------------------------------- linear algebra over F_l
def _rref(A, l):
    """In-place RREF of A (int64, mod l). Returns pivot list [(row, col)]."""
    rows, cols = A.shape
    piv = []
    r = 0
    for c in range(cols):
        if r == rows:
            break
        nz = np.nonzero(A[r:, c])[0]
        if nz.size == 0:
            continue
        p = r + nz[0]
        if p != r:
            A[[r, p]] = A[[p, r]]
        A[r] = (A[r] * pow(int(A[r, c]), l - 2, l)) % l
        col = A[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if nzr.size:
            A[nzr] = (A[nzr] - np.outer(col[nzr], A[r])) % l
        piv.append((r, c)); r += 1
    return piv

def kernel_vector(cols, l):
    """cols: list of m vectors in F_l^rho (np arrays). Nonzero c in F_l^m with sum c_j col_j = 0, or None."""
    m = len(cols)
    if m == 0:
        return None
    A = np.array(cols, dtype=np.int64).T % l
    if A.ndim == 1:
        A = A.reshape(0, m)
    piv = _rref(A, l)
    pc = {c for _, c in piv}
    free = [c for c in range(m) if c not in pc]
    if not free:
        return None
    f = free[0]                      # first free column: support <= rank + 1 and among the earliest columns
    c = np.zeros(m, dtype=np.int64); c[f] = 1
    for r, pcol in piv:
        c[pcol] = (-A[r, f]) % l
    return c

def affine_solve(cols, target, l):
    """c in F_l^m with sum c_j col_j = target, or None."""
    m = len(cols)
    A = np.array(list(cols) + [target], dtype=np.int64).T % l
    piv = _rref(A, l)
    if any(c == m for _, c in piv):
        return None
    c = np.zeros(m, dtype=np.int64)
    for r, pcol in piv:
        c[pcol] = A[r, m]
    return c

# ---------------------------------------------------------------- group bookkeeping
class Group:
    def __init__(self, moduli):
        self.m = [int(x) for x in moduli]
        self.comps = []                               # (index j, prime l, exponent a) with l^a || m_j
        for j, mj in enumerate(self.m):
            x, l = mj, 2
            while x > 1:
                if l * l > x:
                    self.comps.append((j, x, 1)); break
                if x % l == 0:
                    a = 0
                    while x % l == 0:
                        x //= l; a += 1
                    self.comps.append((j, l, a))
                l += 1
        self.primes = sorted({l for _, l, _ in self.comps})

    def add(self, x, y):
        return tuple((a + b) % m for a, b, m in zip(x, y, self.m))

    def zero(self):
        return tuple(0 for _ in self.m)

    def layers(self, order=None):
        """List of (l, k, active comps) for the filtration; order: list of primes."""
        out = []
        for l in (order or self.primes):
            cs = [(j, a) for j, ll, a in self.comps if ll == l]
            for k in range(1, max((a for _, a in cs), default=0) + 1):
                act = [(j, a) for j, a in cs if a >= k]
                out.append((l, k, act))
        return out

def layer_image(G, x, l, k, act):
    v = []
    for j, a in act:
        r = x[j] % (l ** a)
        assert r % (l ** (k - 1)) == 0, "element not in the layer subgroup"
        v.append((r // l ** (k - 1)) % l)
    return np.array(v, dtype=np.int64)

# ---------------------------------------------------------------- elements and streams
class Elt:
    """reals[c] = (sorted tuple of base indices, full value) realising c * (layer value), c = 1..mult."""
    __slots__ = ("v", "reals")
    def __init__(self, v, reals):
        self.v, self.reals = v, reals

def union(G, rs):
    idx, val = [], G.zero()
    for ix, vv in rs:
        idx.extend(ix); val = G.add(val, vv)
    return tuple(sorted(idx)), val

class Stream:
    """Pull interface with push-back buffer."""
    def __init__(self, gen):
        self.gen, self.buf = gen, []
    def pull(self):
        if self.buf:
            return self.buf.pop()
        return self.gen()
    def push(self, e):
        self.buf.append(e)

class Exhausted(Exception):
    pass

def _take(stream, n):
    out = []
    for _ in range(n):
        e = stream.pull()
        if e is None:
            for f in out:
                stream.push(f)
            raise Exhausted
        out.append(e)
    return out

SQRT = {}
def sq_table(l):
    """t -> (a, b) squares with a - nu*b = t, and nu (least nonsquare)."""
    if l in SQRT:
        return SQRT[l]
    sq = sorted({x * x % l for x in range(l)})
    nu = next(x for x in range(2, l) if x not in sq) if l > 2 else None
    tab = {}
    for a in sq:
        for b in sq:
            tab.setdefault((a - nu * b) % l, (a, b))
    for t in range(1, l):                       # prefer single-coset decompositions (smaller supports)
        if t in sq:
            tab[t] = (t, 0)
        else:
            b = (-t * pow(nu, l - 2, l)) % l
            if b in sq:
                tab[t] = (0, b)
    tab[0] = (0, 0)
    SQRT[l] = (nu, sq, tab)
    return SQRT[l]

# ---------------------------------------------------------------- layer stages
class QRRel:
    """Produces relations {elt: coef} with coef in the nonzero squares of F_l, sum coef*v = 0 (layer).
    Norm trick: blocks of rho+1 inputs -> full block elements -> LA over rho+1 of them."""
    def __init__(self, src, l, rho):
        self.src, self.l, self.rho = src, l, rho
        self.nu, self.sq, self.tab = sq_table(l)
        self.blocks = []

    def _block(self):
        l = self.l
        E = _take(self.src, self.rho + 1)
        t = kernel_vector([e.v for e in E], l)
        a, b, used = {}, {}, set()
        for e, ti in zip(E, t):
            ti = int(ti)
            if ti == 0:
                self.src.push(e); continue
            ai, bi = self.tab[ti]
            if ai: a[e] = ai
            if bi: b[e] = bi
        w = np.zeros(self.rho, dtype=np.int64)
        for e, c in b.items():
            w = (w + c * e.v) % l
        return (w, a, b)

    def pull(self):
        l = self.l
        while True:
            # zero-valued block: a ready relation
            for i, (w, a, b) in enumerate(self.blocks):
                if not w.any():
                    del self.blocks[i]
                    return b if b else a
            need = self.rho + 1 - len(self.blocks)
            if need <= 0:
                break
            self.blocks.append(self._block())
        B = self.blocks[: self.rho + 1]
        lam = kernel_vector([w for w, _, _ in B], l)
        rel, keep = {}, []
        for (w, a, b), lk in zip(B, lam):
            lk = int(lk)
            if lk == 0:
                keep.append((w, a, b)); continue
            if lk in self.sq:                                    # lk * w realised by lk * b
                for e, c in b.items(): rel[e] = (lk * c) % l
            else:                                               # lk * w = (lk/nu) * (nu w), realised by (lk/nu) * a
                s = (lk * pow(self.nu, l - 2, l)) % l
                for e, c in a.items(): rel[e] = (s * c) % l
        self.blocks = keep + self.blocks[self.rho + 1:]
        assert rel
        return rel

class Doubler:
    """Ternary relations (coef +-1) over m-fold inputs -> 2m-fold elements (l = 5)."""
    def __init__(self, G, rel_src, l):
        self.G, self.rel, self.l = G, rel_src, l
        self.zeros = []
    def pull(self):
        rel = self.rel.pull()
        Kp = [e for e, c in rel.items() if c == 1]
        Km = [e for e, c in rel.items() if c == self.l - 1]
        assert len(Kp) + len(Km) == len(rel)
        m = min(min(len(e.reals) for e in rel), self.l - 1)
        w = np.zeros(len(next(iter(rel)).v), dtype=np.int64)
        for e in Kp:
            w = (w + e.v) % self.l
        copies = []
        for K in (Kp, Km):
            if not K: continue
            for j in range(m):
                copies.append(union(self.G, [e.reals[j + 1] if j == 0 else _copy(self.G, e, j) for e in K]))
        return _from_copies(self.G, w, copies)

def _copy(G, e, j):
    """j-th disjoint copy of a k-fold element: reals[j+1] minus reals[j] as index sets."""
    hi, lo = e.reals[j + 1], e.reals[j]
    s = set(lo[0]); ix = tuple(i for i in hi[0] if i not in s)
    val = tuple((a - b) % m for a, b, m in zip(hi[1], lo[1], G.m))
    return ix, val

def _from_copies(G, w, copies):
    reals, acc = {}, []
    for c in range(1, len(copies) + 1):
        acc.append(copies[c - 1]); reals[c] = union(G, acc)
    return Elt(w, reals)

class FinalLA:
    """Inputs with mult >= l-1 (all coefficients realisable). One output = subset with zero layer image."""
    def __init__(self, src, l, rho):
        self.src, self.l, self.rho = src, l, rho
    def pull(self, G, forced=None):
        l = self.l
        if forced is None:
            E = []
            while len(E) < self.rho + 1:
                e = self.src.pull()
                if e is None:
                    for f in E: self.src.push(f)
                    raise Exhausted
                if not e.v.any():
                    for f in E: self.src.push(f)
                    return e.reals[1]
                E.append(e)
            lam = kernel_vector([e.v for e in E], l)
        else:
            E = []
            target = (-forced[2]) % l
            lam = None
            if not target.any():
                return forced[0], forced[1]
            while lam is None:
                e = self.src.pull()
                if e is None:
                    for f in E: self.src.push(f)
                    raise Exhausted
                E.append(e)
                if len(E) >= self.rho:
                    lam = affine_solve([f.v for f in E], target, l)
        rs = [] if forced is None else [(forced[0], forced[1])]
        for e, lk in zip(E, lam):
            lk = int(lk)
            if lk:
                rs.append(e.reals[lk])
            else:
                self.src.push(e)
        return union(G, rs)

class L3Final:
    """l = 3: QR relations are already 0/1 (coef 1); targets use the full block elements."""
    def __init__(self, G, qr):
        self.G, self.qr = G, qr
    def pull(self, G, forced=None):
        if forced is None:
            rel = self.qr.pull()
            assert all(c == 1 for c in rel.values())
            return union(G, [e.reals[1] for e in rel])
        target = (-forced[2]) % 3
        if not target.any():
            return forced[0], forced[1]
        qr = self.qr
        while len(qr.blocks) < qr.rho:
            qr.blocks.append(qr._block())
        lam = None
        while lam is None:
            lam = affine_solve([w for w, _, _ in qr.blocks], target, 3)
            if lam is None:
                qr.blocks.append(qr._block())
        rs, keep = [(forced[0], forced[1])], []
        for (w, a, b), lk in zip(qr.blocks, lam):
            lk = int(lk)
            if lk == 0: keep.append((w, a, b)); continue
            mp = b if lk == 1 else a                              # 1*w by b, 2*w = nu*w by a (nu = 2)
            rs.extend(e.reals[1] for e in mp)
        qr.blocks = keep
        return union(G, rs)

# ---------------------------------------------------------------- the layered solver
class Solver:
    def __init__(self, G, X, order=None):
        self.G, self.X = G, [tuple(int(a) % m for a, m in zip(x, G.m)) for x in X]
        self.L = G.layers(order)
        self.pulled = 0
        self._next = 0
        self.levels = []                     # level d stream: elements zero in layers < d (as (idx, val))
        base = Stream(self._base)
        self.levels.append(base)
        for d, (l, k, act) in enumerate(self.L):
            self.levels.append(self._make_level(d, l, k, act))

    def _base(self):
        if self._next >= len(self.X):
            return None
        i = self._next; self._next += 1; self.pulled = self._next
        return (( i,), self.X[i])

    def _make_level(self, d, l, k, act):
        G, rho = self.G, len(act)
        below = self.levels[d]
        def to_elt():
            r = below.pull()
            if r is None:
                return None
            return Elt(layer_image(G, r[1], l, k, act), {1: r})
        inp = Stream(to_elt)
        if l == 2:
            top = FinalLA(inp, 2, rho)
        elif l == 3:
            top = L3Final(G, QRRel(inp, 3, rho))
        elif l == 5:
            dbl1 = Stream(Doubler(G, QRRel(inp, 5, rho), 5).pull)
            dbl2 = Stream(Doubler(G, QRRel(dbl1, 5, rho), 5).pull)
            top = FinalLA(dbl2, 5, rho)
        else:
            raise NotImplementedError(f"layer prime {l}: no polynomial 0/1 method implemented")
        self._tops = getattr(self, "_tops", []) + [(top, l, k, act)]
        def gen():
            try:
                return top.pull(G)
            except Exhausted:
                return None
        return Stream(gen)

    def zero_sum(self):
        """Nonempty I with sum = 0, or None if the pool runs out."""
        r = self.levels[-1].pull()
        if r is None:
            return None
        assert r[1] == self.G.zero()
        return list(r[0])

    def target(self, u):
        """I (possibly empty) with sum_I x + u = 0; u is a forced element carried through all layers."""
        cur = ((), tuple(int(a) % m for a, m in zip(u, self.G.m)))
        for d, (top, l, k, act) in enumerate(self._tops):
            v = layer_image(self.G, cur[1], l, k, act)
            try:
                cur = top.pull(self.G, forced=(cur[0], cur[1], v))
            except Exhausted:
                return None
        assert cur[1] == self.G.zero()
        return list(cur[0])

def bound(G, order=None):
    """Deterministic element count prod f_l(rho) for one zero-sum."""
    f = {2: lambda r: r + 1, 3: lambda r: (r + 1) ** 2, 5: lambda r: (r + 1) ** 5}
    out = 1
    for l, k, act in G.layers(order):
        out *= f[l](len(act)) if l in f else float('inf')
    return out
