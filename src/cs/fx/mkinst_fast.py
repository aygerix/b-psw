"""Same output as e1/mkinst.py (bbc instance files + index.tsv), but reads the record files once.
Usage: T0=... python3 mkinst_fast.py OUTDIR D1[,D2,...] recfile1 [recfile2 ...]"""
import sys, os
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lagr_bb import Inst
from mkinst import write
out = sys.argv[1]; Ds = [int(x) for x in sys.argv[2].split(',')]; files = sys.argv[3:]
os.makedirs(out, exist_ok=True)
pools = defaultdict(dict); bad = defaultdict(int)
for fn in files:
    for l in open(fn):
        a = l.split()
        if len(a) != 7: continue
        D, p, t, c1, c2, o1, o2 = map(int, a)
        key = (D, c1, c2)
        if p in pools[key] and pools[key][p] != (o1, o2): bad[key] += 1
        pools[key][p] = (o1, o2)
idx = open(os.path.join(out, 'index.tsv'), 'a')
for D in Ds:
    for key in sorted((k for k in pools if k[0] == D), key=lambda k: -len(pools[k])):
        inst = Inst(pools[key]); fn = os.path.join(out, f"inst_D{D}_{key[1]}_{key[2]}.txt")
        write(inst, fn)
        idx.write(f"{D}\t{key[1]}\t{key[2]}\t{inst.n}\t{len(inst.items)}\t{len(inst.rows)}\t{bad[key]}\t{fn}\n")
idx.close()
