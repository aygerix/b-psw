"""Exact M(lambda) by exhaustive enumeration of the sides of the conflict rows (rows needed on both sides by some prime),
each pattern solved as a max-weight closure with scipy's maximum_flow (independent of bbc/bb2: no branch-and-bound, no
deletion lemma).  Runs on every certificate instance with <= CMAX conflict rows; compares with the claimed M_ub.
Usage: python3 bfx.py CMAX runfiles..."""
import sys, os, math, itertools
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_flow
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from instlib import read_inst, item_q
CMAX = int(sys.argv[1]); runs = sys.argv[2:]
claimed = {}
for fn in runs:
    cur = None
    for l in open(fn):
        a = l.split()
        if not a: continue
        if a[0] == '==': cur = (int(a[1]), int(a[2]), int(a[3]))
        elif a[0] == 'lam' and len(a) > 4 and a[2] == 'DONE' and cur: claimed.setdefault(cur, {})[round(float(a[1]), 6)] = float(a[4])
B = os.path.expanduser('~/b-psw/src/cs/')
insts = []
for idx, root in ((B + 'e1/out_index_copy.tsv', B + 'e1/'), (B + 'fx/instxD/index.tsv', B + 'fx/')):
    for l in open(idx):
        a = l.split('\t'); insts.append(((int(a[0]), int(a[1]), int(a[2])), root + a[7].strip()))
KS = 1 << 16            # int32-safe (scipy maximum_flow works in int32): total source capacity <= 6523 * 2^16 < 2^31
from scipy.sparse.csgraph import breadth_first_order
def closure_value(alive, prim, items, cost, lam):
    """max-weight closure with costs rounded down to 1/KS (relaxed value >= true); returns (rounded bound, exact value of
    the minimal optimal closure found = a feasible lower bound)."""
    ia = [i for i in range(len(prim)) if alive[i]]
    if not ia: return 0.0, 0.0
    used = set()
    for i in ia:
        for j in prim[i][1]:
            while j >= 0 and j not in used: used.add(j); j = items[j][3]
    used = sorted(used); inode = {j: 2 + len(ia) + k for k, j in enumerate(used)}
    r, c, w = [], [], []
    BIG = 1 << 30
    for k, i in enumerate(ia):
        r.append(0); c.append(2 + k); w.append(KS)
        for j in prim[i][1]: r.append(2 + k); c.append(inode[j]); w.append(BIG)
    for j in used:
        if items[j][3] >= 0: r.append(inode[j]); c.append(inode[items[j][3]]); w.append(BIG)
        cc = int(math.floor(lam * cost[j] * KS))
        if cc > 0: r.append(inode[j]); c.append(1); w.append(cc)
    n = 2 + len(ia) + len(used)
    G = csr_matrix((np.array(w, dtype=np.int32), (np.array(r), np.array(c))), shape=(n, n))
    res = maximum_flow(G, 0, 1)
    ub = (len(ia) * KS - res.flow_value) / KS
    R = (G - res.flow).tocsr(); R.data[R.data < 0] = 0; R.eliminate_zeros()
    reach = set(breadth_first_order(R, 0, directed=True, return_predecessors=False).tolist())
    selp = [ia[k] for k in range(len(ia)) if 2 + k in reach]
    u = set()
    for i in selp:
        for j in prim[i][1]:
            while j >= 0 and j not in u: u.add(j); j = items[j][3]
    return ub, len(selp) - lam * sum(cost[j] for j in u)
done = 0; bad = []; skipped = 0; maxdiff = 0
for key, fn in insts:
    nP, nI, nR, items, prim = read_inst(fn); qs = item_q(items)
    cost = [math.log2(qs[j] - 1) if items[j][2] == 1 else math.log2(qs[j]) for j in range(nI)]
    sides = {}
    for p, its in prim:
        for j in its: sides.setdefault(items[j][0], set()).add(items[j][1])
    conf = sorted(r for r, s in sides.items() if len(s) == 2)
    if len(conf) > CMAX: skipped += 1; continue
    need = [{(items[j][0], items[j][1]) for j in its} for p, its in prim]
    for lam in (0.3, 0.25, 0.2):
        best = 0.0; bestlb = 0.0
        for pat in itertools.product((0, 1), repeat=len(conf)):
            fix = dict(zip(conf, pat))
            alive = [all(not (r in fix and fix[r] != s) for (r, s) in nd) for nd in need]
            ub_, lb_ = closure_value(alive, prim, items, cost, lam)
            best = max(best, ub_); bestlb = max(bestlb, lb_)
        cl = claimed.get(key, {}).get(lam)
        done += 1
        d = (cl - best) if cl is not None else None
        # exact M lies in [bestlb, best]; best - bestlb <= (#items)/KS.  Claimed M_ub must be >= bestlb (else claim is too
        # small = WRONG) and should be <= best + 1e-5 (else claim is loose).
        if cl is None or bestlb > cl + 1e-6 or cl > best + 1e-5:
            bad.append((key, lam, len(conf), round(bestlb, 6), round(best, 6), cl))
        else: maxdiff = max(maxdiff, abs(cl - best))
print(f"instances with <= {CMAX} conflict rows: {done // 3} (x 3 lambdas = {done}); skipped {skipped}; "
      f"disagreements (exact > claimed, or claimed > exact + 1e-5, or missing): {len(bad)}; max |claimed - exact| otherwise {maxdiff:.2e}")
for b in bad[:30]: print("  ", b)
