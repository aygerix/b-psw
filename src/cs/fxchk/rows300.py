"""Count, in the T0=300 instance D=-15 (2,0), how many primes need each row q in (200,300] on each side."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from instlib import read_inst, item_q
nP, nI, nR, items, prim = read_inst(os.path.expanduser('~/b-psw/src/cs/fx/inst300/inst_D-15_2_0.txt')); qs = item_q(items)
cnt = {}
for p, its in prim:
    for j in its:
        q = qs[j] if items[j][2] == 1 else qs[j]
        if q > 200: cnt[(q, items[j][1])] = cnt.get((q, items[j][1]), 0) + 1
rows = sorted({q for q, s in cnt})
print("new rows:", len(rows), [(q, cnt.get((q, 0), 0), cnt.get((q, 1), 0)) for q in rows])
v = [cnt.get((q, s), 0) for q in rows for s in (0, 1)]; print("per side min/max:", min(v), max(v), " per row total min/max:", min(cnt.get((q,0),0)+cnt.get((q,1),0) for q in rows), max(cnt.get((q,0),0)+cnt.get((q,1),0) for q in rows))
