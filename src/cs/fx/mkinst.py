"""Write bbc.c instance files for every (D, 2-adic class) of the E1 pool.
Usage: python3 mkinst.py OUTDIR D1[,D2,...] recfile1 [recfile2 ...]   -> OUTDIR/inst_D<D>_<j1>_<j2>.txt + index.tsv
Instance format: see bbc.c.  Rows are indexed in increasing q; items (row, side, level) with incremental cost
log2(q-1) at level 1 and log2(q) above; each prime lists the top item of each of its requirements."""
import sys, os
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lagr_bb import load, Inst, cost

def write(inst, fn):
    rid = {q: k for k, q in enumerate(inst.rows)}
    items = sorted(inst.items.items(), key=lambda t: t[1])
    with open(fn, 'w') as f:
        f.write(f"{inst.n} {len(items)} {len(inst.rows)}\n")
        for (q, s, e), j in items:
            par = inst.items[(q, s, e - 1)] if e > 1 else -1
            f.write(f"{rid[q]} {s} {e} {par} {cost(q, e):.17g}\n")
        for p, req in zip(inst.P, inst.req):
            f.write(f"{p} {len(req)} " + " ".join(str(inst.items[(q, s, a)]) for q, s, a in req) + "\n")

def main():
    out = sys.argv[1]; Ds = [int(x) for x in sys.argv[2].split(',')]; files = sys.argv[3:]
    os.makedirs(out, exist_ok=True)
    recs = {}
    for fn in files:
        for l in open(fn):
            a = l.split()
            if len(a) == 7: recs.setdefault(int(a[0]), {})[int(a[1])] = (int(a[3]), int(a[4]))
    idx = open(os.path.join(out, 'index.tsv'), 'a')
    for D in Ds:
        for c, cnt in sorted(Counter(recs.get(D, {}).values()).items(), key=lambda t: -t[1]):
            pool, bad = load(D, c[0], c[1], files)
            inst = Inst(pool)
            fn = os.path.join(out, f"inst_D{D}_{c[0]}_{c[1]}.txt")
            write(inst, fn)
            idx.write(f"{D}\t{c[0]}\t{c[1]}\t{inst.n}\t{len(inst.items)}\t{len(inst.rows)}\t{bad}\t{fn}\n")
    idx.close()

if __name__ == "__main__":
    main()
