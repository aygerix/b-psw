"""Split a bbc instance by fixing the sides of given rows: for every side pattern s in {F, L}^k, delete the primes that
need one of these rows on the other side.  Every configuration has some side pattern on these rows and is feasible, with the
same N and bits, in the matching sub-instance, so M(lambda) = max over sub-instances (certificate = max of their M_ub).
Usage: python3 split_inst.py INST OUTPREFIX q1[,q2,...]     -> OUTPREFIX_<sides>.txt   (sides string like FL)"""
import sys, itertools
fn, pre, qs = sys.argv[1], sys.argv[2], [int(x) for x in sys.argv[3].split(',')]
L = open(fn).read().split('\n')
nP, nI, nR = map(int, L[0].split())
items = [l.split() for l in L[1:1 + nI]]
# row index -> q: rows sorted increasing q in mkinst; recover q from the cost of level-1 items (log2(q - 1)) -- check
import math
rowq = {}
for it in items:
    r, sd, lev, par, c = int(it[0]), int(it[1]), int(it[2]), int(it[3]), float(it[4])
    if lev == 1: rowq[r] = round(2 ** c) + 1
    else: assert abs(2 ** c - round(2 ** c)) < 1e-6
qrow = {q: r for r, q in rowq.items()}
assert sorted(rowq.values()) == [rowq[r] for r in sorted(rowq)], "rows not increasing in q"
prim = [l.split() for l in L[1 + nI:1 + nI + nP]]
for pat in itertools.product((0, 1), repeat=len(qs)):
    fix = {qrow[q]: s for q, s in zip(qs, pat)}
    keep = []
    for pl in prim:
        ok = True
        for j in map(int, pl[2:]):
            r, sd = int(items[j][0]), int(items[j][1])
            if r in fix and fix[r] != sd: ok = False; break
        if ok: keep.append(pl)
    out = pre + '_' + ''.join('FL'[s] for s in pat) + '.txt'
    with open(out, 'w') as f:
        f.write(f"{len(keep)} {nI} {nR}\n")
        for it in L[1:1 + nI]: f.write(it + '\n')
        for pl in keep: f.write(' '.join(pl) + '\n')
    print(out, len(keep))
