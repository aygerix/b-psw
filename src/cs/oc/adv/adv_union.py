#!/usr/bin/env python3
"""Adversarial extension (line OC): CROSS-FAMILY union classes. The line searched each family separately; here all
usable members (both odd orders 1500-smooth, values recomputed by adv_rows.py) of all families are pooled per
(D, j1, j2), deduplicated by p, OC3-dead classes (j1, j2 >= 2) dropped, and searched:
  exhaustive if the class has <= EXH members, else all subsets of size <= KMAX.
Necessary conditions (Lemma 1 for this D): prod eps = -1, odd ord_p 2 | n - 1, odd z(p) | n + 1, v2 room.
Survivors -> src/bpsw.py. Control: on subsets with n < 2^CTLB the Fermat half is compared with strong_fermat2(n)
and, when selfridge_D(n) == D, the Lucas half with strong_lucas(n).
Usage: adv_union.py EXH KMAX CTLB rows_file [rows_file ...]"""
import sys, os, time, itertools
from collections import defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))
from bpsw import bpsw_pass, strong_fermat2, strong_lucas, selfridge_D
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from adv_rows import member

def v2(x): return (x & -x).bit_length() - 1
EXH, KMAX, CTLB = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
cls = defaultdict(dict)
for path in sys.argv[4:]:
    for line in open(path):
        t = line.split()
        if len(t) != 13: continue
        fam, par, D, eps, j1, j2, oo, oz = t[0], int(t[1]), int(t[3]), int(t[4]), int(t[5]), int(t[6]), int(t[9]), int(t[10])
        if oo <= 0 or oz <= 0: continue
        p = member(fam, par)
        cls[(D, j1, j2)][p] = (p, oo, oz, eps)
t0 = time.time()
tot_sub = tot_pass = tot_hit = 0
ctl = dict(n=0, fagree=0, fpos=0, lcmp=0, lagree=0, lpos=0)
nclass = 0
for (D, j1, j2), mem in sorted(cls.items(), key=lambda kv: -len(kv[1])):
    if j1 >= 2 and j2 >= 2: continue           # OC3: no inert prime possible
    M = sorted(mem.values())
    m = len(M)
    if m < 2: continue
    nclass += 1
    kmax = m if m <= EXH else KMAX
    sub = pas = 0
    for k in range(2, kmax + 1):
        for S in itertools.combinations(range(m), k):
            sub += 1
            e = 1; n = 1
            for i in S: e *= M[i][3]; n *= M[i][0]
            fh = all((n - 1) % M[i][1] == 0 for i in S) and j1 <= v2(n - 1)
            lh = all((n + 1) % M[i][2] == 0 for i in S) and j2 <= v2(n + 1)
            if n.bit_length() <= CTLB:
                ctl['n'] += 1; sf = strong_fermat2(n); ctl['fagree'] += (sf == fh); ctl['fpos'] += sf
                if sf != fh: print('FERMAT-HALF DISAGREE', D, (j1, j2), [M[i][0] for i in S], fh, sf)
                if selfridge_D(n) == D:
                    sl = strong_lucas(n); ctl['lcmp'] += 1; ctl['lagree'] += (sl == lh); ctl['lpos'] += sl
                    if sl != lh: print('LUCAS-HALF DISAGREE', D, (j1, j2), [M[i][0] for i in S], lh, sl)
            if e != -1 or not fh or not lh: continue
            pas += 1
            ok = bpsw_pass(n)
            print(f'  NECESSARY-CONDITION PASS D={D} class=({j1},{j2}) S={[M[i][0] for i in S]} selD={selfridge_D(n)} bpsw={ok}', flush=True)
            if ok: tot_hit += 1; print('*** BPSW HIT', [M[i][0] for i in S], flush=True)
    tot_sub += sub; tot_pass += pas
    if m >= 6 or pas: print(f'D={D} class=({j1},{j2}) members={m} {"exhaustive" if m <= EXH else "k<=%d" % KMAX} subsets={sub} pass={pas}  {time.time()-t0:.1f}s', flush=True)
print(f'classes searched {nclass}, subsets {tot_sub}, necessary-condition passes {tot_pass}, BPSW hits {tot_hit}, {time.time()-t0:.1f}s')
print('control:', ctl)
