"""Audit of the E1 pool and instance files (independent of lagr_bb.py / mkinst.py code).
1. Union of census3 (p < 1e10) and CHM records: total (D, p) records, conflicting duplicates, classes, admissible classes.
2. For every instance listed in index.tsv: rebuild prime -> {(q, side, a)} from the records and compare with the instance
   file semantics (items -> (row q, side, level); prime line -> top items).
Usage: python3 aud_pool.py c3file chmfile indexfile instdir"""
import sys, os
from collections import defaultdict
PR = [q for q in range(3, 201) if all(q % d for d in range(2, int(q ** .5) + 1))]
def fac(x):
    f = {}
    for q in PR:
        while x % q == 0: x //= q; f[q] = f.get(q, 0) + 1
    return f if x == 1 else None
c3f, chmf, idxf, idir = sys.argv[1:5]
rec, conflicts, nonsmooth, src = {}, 0, 0, defaultdict(int)
for fn, tag in ((c3f, 'c3'), (chmf, 'chm')):
    for l in open(fn):
        t = tuple(map(int, l.split()))
        if len(t) != 7: continue
        D, p, s, j1, j2, o1, o2 = t
        if fac(o1) is None or fac(o2) is None: nonsmooth += 1; continue
        k = (D, p)
        if k in rec and rec[k] != t: conflicts += 1
        if k not in rec: src[tag] += 1
        rec[k] = t
cls = defaultdict(list)
for (D, p), t in rec.items(): cls[(D, t[3], t[4])].append(t)
adm = [c for c in cls if min(c[1], c[2]) <= 1]
inadm_inert = sum(1 for c in cls if min(c[1], c[2]) >= 2 for t in cls[c] if t[2] == 0)
print(f"union records {len(rec)} (new from c3 {src['c3']}, new from CHM {src['chm']}), conflicting duplicates {conflicts}, "
      f"non-smooth lines {nonsmooth}; classes {len(cls)}, admissible (min(j1,j2) <= 1) {len(adm)}, "
      f"inert primes in inadmissible classes {inadm_inert}")
print("records per D:", {D: sum(1 for k in rec if k[0] == D) for D in sorted({k[0] for k in rec})})
# instance comparison
bad_inst = 0; ninst = 0
for l in open(idxf):
    D, j1, j2, n, ni, nr, bad, fn = l.rstrip('\n').split('\t'); D, j1, j2 = int(D), int(j1), int(j2)
    fn = os.path.join(idir, os.path.basename(fn))
    mine = {}
    for t in cls.get((D, j1, j2), []):
        req = frozenset([(q, 0, a) for q, a in fac(t[5]).items()] + [(q, 1, a) for q, a in fac(t[6]).items()])
        mine[t[1]] = req
    f = open(fn); nP, nI, nR = map(int, f.readline().split())
    items = [f.readline().split() for _ in range(nI)]
    rowsq = sorted({q for req in mine.values() for (q, s, a) in req})
    theirs = {}
    for _ in range(nP):
        a = f.readline().split(); p = int(a[0]); k = int(a[1])
        req = set()
        for x in a[2:2 + k]:
            r, s, e, par, c = items[int(x)]
            req.add((rowsq[int(r)], int(s), int(e)))
        theirs[p] = frozenset(req)
    # item cost check
    import math
    for r, s, e, par, c in items:
        q = rowsq[int(r)]; e = int(e)
        cc = math.log2(q - 1) if e == 1 else math.log2(q)
        if abs(cc - float(c)) > 1e-12: print("  cost mismatch", fn, q, e, c, cc); bad_inst += 1; break
    ninst += 1
    if mine != theirs or nP != len(mine) or len(rowsq) != nR:
        bad_inst += 1
        print(f"  MISMATCH {fn}: mine {len(mine)} theirs {len(theirs)} rows {len(rowsq)} vs {nR}; "
              f"diff {len(set(mine.items()) ^ set(theirs.items()))}")
print(f"instances compared {ninst}, mismatching {bad_inst}")
print("instances in index but not in my classes:", ninst - len(cls) if ninst != len(cls) else 0, "| my classes", len(cls))
