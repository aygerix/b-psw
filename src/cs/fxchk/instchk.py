"""Verify every certificate instance file against the raw records (own builder).  For each (D, class) instance listed in
the index files: same prime set as the record union restricted to (D, j1, j2) (T0-smooth orders), same requirement
multiset per prime, item-chain sanity.  Also verify the split sub-instances (split15, split11) cover the full instance."""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from instlib import *
B = os.path.expanduser('~/b-psw/src/cs/')
def load(files):
    pool = {}
    for fn in files:
        with open(fn) as f:
            for l in f:
                a = l.split()
                if len(a) == 7:
                    r = tuple(int(x) for x in a); pool[(r[0], r[1])] = r
    return pool
pool = load([B + 'e1/c3_1e10.txt', B + 'e1/chm200_pool.txt', B + 'fx/c3_xD_1e10.txt', B + 'fx/chm200_pool_xD.txt'])
byc = {}
for (D, p), r in pool.items(): byc.setdefault((D, r[3], r[4]), {})[p] = (r[5], r[6])
def verify(fn, D, j1, j2, T0=200, quiet=True):
    nP, nI, nR, items, prim = read_inst(fn)
    reqs, qs = inst_reqs(items, prim)
    bad = check_items(items, qs)
    mine = {}
    nonsmooth = 0
    for p, (o1, o2) in byc.get((D, j1, j2), {}).items():
        rq = req_from_record(o1, o2, T0)
        if rq is None: nonsmooth += 1; continue
        mine[p] = rq
    missing = set(mine) - set(reqs); extra = set(reqs) - set(mine)
    diff = sum(1 for p in set(mine) & set(reqs) if mine[p] != reqs[p])
    return dict(n=nP, mine=len(mine), missing=len(missing), extra=len(extra), diff=diff, itembad=bad, nonsmooth=nonsmooth), reqs, items, qs
tot = 0; badl = []
for idxfn, root in ((B + 'e1/out_index_copy.tsv', B + 'e1/'), (B + 'fx/instxD/index.tsv', B + 'fx/')):
    if not os.path.exists(idxfn): print("missing index", idxfn); continue
    for l in open(idxfn):
        a = l.split('\t')
        D, j1, j2, n = int(a[0]), int(a[1]), int(a[2]), int(a[3]); fn = root + a[7].strip()
        res, *_ = verify(fn, D, j1, j2)
        tot += 1
        if res['missing'] or res['extra'] or res['diff'] or res['itembad'] or res['n'] != n or res['mine'] != n: badl.append((D, j1, j2, res))
print(f"instances verified: {tot}; with any mismatch: {len(badl)}")
for b in badl[:20]: print("  MISMATCH", b)
# classes present in records but with no instance file
listed = set()
for idxfn in (B + 'e1/out_index_copy.tsv', B + 'fx/instxD/index.tsv'):
    for l in open(idxfn):
        a = l.split('\t'); listed.add((int(a[0]), int(a[1]), int(a[2])))
print("record classes without an instance:", sorted(set(byc) - listed)[:20], len(set(byc) - listed))
# split coverage
def split_check(full, subs, D, j1, j2):
    res, reqs, items, qs = verify(full, D, j1, j2)
    full_items = open(full).read().split('\n')[1:1 + res['n'] and None]
    rows = sorted(set(q for p in reqs for (q, s, a) in reqs[p]))
    need = {p: {(q, s) for (q, s, a) in reqs[p]} for p in reqs}
    out = []
    for name, sfn in subs:
        sn, sI, sR, sitems, sprim = read_inst(sfn)
        same_items = (sitems == items)
        sp = {p for p, _ in sprim}
        # find the side pattern: fixed (q, s) such that sp == {p : no need (q, 1-s)} ; search over <= 4 rows among small q
        found = None
        cand = rows[:8]
        for k in range(1, 5):
            for rs in itertools.combinations(cand, k):
                for pat in itertools.product((0, 1), repeat=k):
                    comp = {p for p in reqs if all((q, 1 - s) not in need[p] for q, s in zip(rs, pat))}
                    if comp == sp: found = (rs, ''.join('FL'[s] for s in pat)); break
                if found: break
            if found: break
        out.append((name, sn, same_items, found))
    return res, out
res, out = split_check(B + 'e1/inst/inst_D-15_2_0.txt', [(s, B + f'fx/split15/m15_{s}.txt') for s in ['FF', 'FL', 'LF_F', 'LF_L', 'LL_F_F', 'LL_F_L', 'LL_L']], -15, 2, 0)
print("D=-15 (2,0) full:", res); [print("   sub", o) for o in out]
res, out = split_check(B + 'e1/inst/inst_D-11_2_0.txt', [(s, B + f'fx/split11/m11_20_{s}.txt') for s in ['F', 'L']], -11, 2, 0)
print("D=-11 (2,0) full:", res); [print("   sub", o) for o in out]
