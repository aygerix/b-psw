"""E1 certificate: exact concave majorant of the small-b frontier N_max(b), per Selfridge D and 2-adic class.

Model (rows <= T0).  A configuration C = (side s_q in {F, L} and cap e_q >= 0 for every odd prime q <= T0).
Its group G(C) = prod_q (Z/q^e_q)^*, bits(C) = sum_q [e_q > 0] (log2(q-1) + (e_q - 1) log2 q).
Pool prime p (record "D p t j1 j2 o1 o2", class (j1, j2)) is USABLE in C iff every q^a || o1 has s_q = F, e_q >= a,
and every q^a || o2 has s_q = L, e_q >= a.  N(C) = number of usable pool primes.
N_max(b) = max {N(C) : bits(C) <= b}  (the ILP of cs_rank.md E1: x_p <= y_{q,side,e}, one side per row, bits <= b).

Certificate.  For lambda > 0 let M(lambda) = max_C [N(C) - lambda bits(C)].  Then for every b,
      N_max(b) <= lambda b + M(lambda)            (weak duality, valid for every lambda),
and if C* attains M(lambda) then N_max(bits(C*)) = N(C*) exactly (C* is on the upper concave hull).
For FIXED sides, max_C [N - lambda bits] is a maximum-weight closure (prime -> top item (q, s, a), item (q, s, e) ->
(q, s, e-1); weights +1 and -lambda*cost), solved exactly by one min cut (its LP is integral).  With some rows' sides
free, giving free rows BOTH side-chains is a relaxation (upper bound).  Branch-and-bound over sides: a node is closed
when its relaxed optimum uses no free row on both sides (then it is a feasible configuration: exact), or when its
bound <= incumbent.  Capacities are integers: prime weight K, item cost floor(lambda c K) rounded DOWN, so every
computed bound is >= the true relaxed value (conservative); incumbents are re-evaluated in exact floating point.
Usage: python3 lagr_bb.py D j1 j2 lam1[,lam2,...] recfile1 [recfile2 ...]   [env T0=200, PMAX=inf, NODECAP=1e6]
"""
import sys, os, math, time
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_flow, breadth_first_order

T0 = int(os.environ.get('T0', 200))
PMAX = float(os.environ.get('PMAX', 'inf'))
NODECAP = int(float(os.environ.get('NODECAP', 1e6)))
K = 1 << 16
INF = (1 << 31) - 1
PR = [q for q in range(3, T0 + 1) if all(q % d for d in range(2, int(q ** .5) + 1))]

def fac(x):
    f = {}
    for q in PR:
        if x % q == 0:
            e = 0
            while x % q == 0: x //= q; e += 1
            f[q] = e
    assert x == 1, "order not T0-smooth"
    return f

def load(D, j1, j2, files):
    """union of records for (D, class); returns dict p -> (o1, o2); reports conflicting duplicates."""
    pool, bad = {}, 0
    for fn in files:
        for l in open(fn):
            a = l.split()
            if len(a) != 7 or int(a[0]) != D: continue
            p, t, c1, c2, o1, o2 = map(int, a[1:])
            if (c1, c2) != (j1, j2) or p > PMAX: continue
            if p in pool and pool[p] != (o1, o2): bad += 1
            pool[p] = (o1, o2)
    return pool, bad

def cost(q, e): return math.log2(q - 1) if e == 1 else math.log2(q)

class Inst:
    def __init__(self, pool):
        self.P = sorted(pool)
        self.req = []                     # per prime: list of (q, side, a), side 0 = F, 1 = L
        emax = {}
        for p in self.P:
            o1, o2 = pool[p]
            r = [(q, 0, a) for q, a in fac(o1).items()] + [(q, 1, a) for q, a in fac(o2).items()]
            self.req.append(r)
            for q, s, a in r: emax[(q, s)] = max(emax.get((q, s), 0), a)
        self.items = {}                   # (q, s, e) -> index
        for (q, s), E in sorted(emax.items()):
            for e in range(1, E + 1): self.items[(q, s, e)] = len(self.items)
        self.rows = sorted({q for q, s in emax})
        self.n = len(self.P)
        self.need = {}                    # (q, s) -> array of prime indices needing row q on side s
        for i, r in enumerate(self.req):
            for q, s, a in r: self.need.setdefault((q, s), []).append(i)
        self.need = {k: np.array(v) for k, v in self.need.items()}
        self.bits_exact = lambda used: sum(cost(q, e) for (q, s, e) in used)

    def closure(self, alive, lam):
        """max-weight closure on primes with alive[i]; returns (ub, selected prime idx, used items set)."""
        idx = np.nonzero(alive)[0]
        if len(idx) == 0: return 0.0, [], set(), [], set()
        nP = len(idx); it = self.items; nI = len(it)
        S, T = 0, 1
        rows_, cols_, caps_ = [], [], []
        pnode = {i: 2 + k for k, i in enumerate(idx)}
        base = 2 + nP
        for k, i in enumerate(idx):
            rows_.append(S); cols_.append(2 + k); caps_.append(K)
            for q, s, a in self.req[i]:
                rows_.append(2 + k); cols_.append(base + it[(q, s, a)]); caps_.append(INF)
        for (q, s, e), j in it.items():
            if e > 1: rows_.append(base + j); cols_.append(base + it[(q, s, e - 1)]); caps_.append(INF)
            c = int(math.floor(lam * cost(q, e) * K * (1 - 1e-12)))
            rows_.append(base + j); cols_.append(T); caps_.append(c)
        N = base + nI
        C = csr_matrix((np.array(caps_, dtype=np.int32), (np.array(rows_), np.array(cols_))), shape=(N, N))
        res = maximum_flow(C, S, T, method='dinic')
        ub = (nP * K - res.flow_value) / K
        Rm = (C - res.flow).tocsr(); Rm.data[Rm.data < 0] = 0; Rm.eliminate_zeros()
        reach = breadth_first_order(Rm, S, directed=True, return_predecessors=False)
        rs = set(reach.tolist())                       # MINIMAL optimal closure (source side reachable from S)
        sel = [i for k, i in enumerate(idx) if 2 + k in rs]
        inv = {j: key for key, j in it.items()}
        used = {inv[n - base] for n in rs if n >= base}
        # MAXIMAL optimal closure = nodes that cannot reach T in the residual graph.  By monotonicity of maximal
        # min-cut source sets under decreasing source capacities (deleting primes), a prime outside it belongs to no
        # optimal closure of any descendant node, so it can be deleted from the whole subtree.
        toT = breadth_first_order(Rm.T.tocsr(), T, directed=True, return_predecessors=False)
        tt = set(toT.tolist())
        selmax = [i for k, i in enumerate(idx) if 2 + k not in tt]
        usedmax = {inv[n - base] for n in range(base, N) if n not in tt and (n - base) in inv}
        return ub, sel, used, selmax, usedmax

    def exact_value(self, sel, lam):
        used = set()
        for i in sel:
            for q, s, a in self.req[i]:
                for e in range(1, a + 1): used.add((q, s, e))
        b = sum(cost(q, e) for q, s, e in used)
        return len(sel) - lam * b, len(sel), b, used

def bb(inst, lam, log=None):
    """returns dict: M_ub (certified), incumbent (value, N, bits, sel), nodes, exact-closed, pruned."""
    t0 = time.time()
    M_ub = -1e300
    stats = dict(nodes=0, exact=0, pruned=0)
    best = (0.0, 0, 0.0, [])               # empty configuration
    stack = [(dict(), np.ones(inst.n, dtype=bool))]   # (fixed sides {q: s}, alive primes)
    while stack:
        fixed, alive = stack.pop()
        stats['nodes'] += 1
        if stats['nodes'] > NODECAP:
            return dict(status='NODECAP', **stats, M_ub=None, best=best, time=time.time() - t0)
        ub, sel, used, selmax, usedmax = inst.closure(alive, lam)
        if ub <= best[0] + 1e-9:
            stats['pruned'] += 1; M_ub = max(M_ub, ub); continue
        def conflicts(us):
            both = {}
            for (q, s, e) in us:
                if e == 1 and q not in fixed: both.setdefault(q, set()).add(s)
            return [q for q, ss in both.items() if len(ss) == 2]
        if not conflicts(used):
            v, n, b, _ = inst.exact_value(sel, lam)
            stats['exact'] += 1; M_ub = max(M_ub, ub)
            if v > best[0]: best = (v, n, b, sel)
            continue
        conflict = conflicts(usedmax)
        alive2 = np.zeros(inst.n, dtype=bool); alive2[selmax] = True    # subtree-safe deletion
        selset = alive2
        def w(q, s):
            k = (q, s)
            return int(selset[inst.need[k]].sum()) if k in inst.need else 0
        q = max(conflict, key=lambda q: (min(w(q, 0), w(q, 1)), q))
        sF, sL = w(q, 0), w(q, 1)
        first = 0 if sF >= sL else 1
        for s in (1 - first, first):        # stack: push weaker first, explore stronger first
            f2 = dict(fixed); f2[q] = s
            a2 = alive2.copy()
            k = (q, 1 - s)
            if k in inst.need: a2[inst.need[k]] = False
            stack.append((f2, a2))
    M_ub = max(M_ub, best[0])
    return dict(status='DONE', **stats, M_ub=M_ub, best=best, time=time.time() - t0)

def main():
    D, j1, j2 = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    lams = [float(x) for x in sys.argv[4].split(',')]
    pool, bad = load(D, j1, j2, sys.argv[5:])
    inst = Inst(pool)
    print(f"D = {D}, class ({j1},{j2}), T0 = {T0}, PMAX = {PMAX:g}: pool {inst.n} primes "
          f"(conflicting duplicate records: {bad}), rows {len(inst.rows)}, items {len(inst.items)}", flush=True)
    for lam in lams:
        r = bb(inst, lam)
        v, n, b, sel = r['best']
        print(f"lambda = {lam:.4f}: {r['status']} M_ub = {r['M_ub'] if r['M_ub'] is None else round(r['M_ub'], 5)} "
              f"best = {v:.5f} at (bits {b:.2f}, N {n}); nodes {r['nodes']} (exact {r['exact']}, pruned {r['pruned']}); "
              f"{r['time']:.1f} s", flush=True)
        if sel:
            print(f"   hull point primes: {[inst.P[i] for i in sel][:60]}", flush=True)

if __name__ == "__main__":
    main()
