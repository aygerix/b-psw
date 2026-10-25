#!/usr/bin/env python3
"""Adversarial recomputation of oc_summary.py statistics (new code) from a row file set, plus inheritance bands.
Usage: adv_summary.py rows_dir inherit_file
rows_dir holds <fam>.rows (adv_rows.py output, recomputed values)."""
import sys, math, statistics, os
from collections import defaultdict
FAMS = [('mer', 'F'), ('wag', 'F'), ('gm+', 'F'), ('gm-', 'F'), ('cyc', 'F'), ('fib', 'L'), ('luc', 'L'), ('fcy', 'L')]
SRC = {'mer': 'mer', 'wag': 'wag', 'gm+': 'gm', 'gm-': 'gm', 'cyc': 'cyc', 'fib': 'fib', 'luc': 'luc', 'fcy': 'fcy'}
def odd(x):
    while x % 2 == 0: x //= 2
    return x
def forced_bits(fam, par):
    if fam in ('cyc', 'fcy'): o = odd(par); return math.log2(o) if o > 1 else 0.0
    return math.log2(par)
d = sys.argv[1]
for D in (5, -7):
    print(f'=== D={D}: fam n inert #cls maxcls maxcls_u slope cost_med cost_min')
    for fam, forced in FAMS:
        R = []
        for line in open(os.path.join(d, SRC[fam] + '.rows')):
            t = line.split()
            if t[0] != fam or int(t[3]) != D: continue
            R.append(dict(par=int(t[1]), bits=int(t[2]), eps=int(t[4]), j=(int(t[5]), int(t[6])), oo=int(t[9]), oz=int(t[10]), rF=float(t[11]), rL=float(t[12])))
        if not R: continue
        cl = defaultdict(list)
        for r in R: cl[r['j']].append(r)
        maxc = max(len(v) for v in cl.values()); maxu = max(sum(1 for r in v if r['oo'] > 0 and r['oz'] > 0) for v in cl.values())
        free = 'rL' if forced == 'F' else 'rF'
        sel = [r for r in R if r['bits'] >= 64 and (forced == 'L' or r['eps'] == -1)]
        slope = sum(r[free] for r in sel) / sum(r['bits'] for r in sel) if sel else float('nan')
        costs = []
        for r in sel:
            if r['oo'] > 0 and r['oz'] > 0: costs.append(math.log2(r['oo']) + math.log2(r['oz']))
            else: costs.append(forced_bits(fam, r['par']) + r[free] + (r['rL'] if (forced == 'L' and D != 5) else 0.0))
        cm = statistics.median(costs) if costs else float('nan'); cmin = min(costs) if costs else float('nan')
        print(f'{fam:4s} {len(R):3d} {sum(1 for r in R if r["eps"] == -1):4d} {len(cl):4d} {maxc:4d} {maxu:4d} {slope:6.3f} {cm:7.1f} {cmin:7.1f}')
# inheritance
rows = [l.split() for l in open(sys.argv[2]) if l.startswith('INH')]
print('=== inheritance rows', len(rows))
for fam in ('cyc', 'fcy'):
    R = [(int(t[3]), float(t[4]), float(t[5]), float(t[6]), float(t[7])) for t in rows if t[1] == fam]
    print(f' {fam}: n={len(R)}; band: n, mean[sm(p+1)-ctl], mean[sm(p-1)-ctl], mean bits')
    for lo, hi in ((0, 64), (64, 128), (128, 256), (256, 512), (512, 1024), (1024, 4096)):
        B = [r for r in R if lo <= r[0] < hi]
        if not B: continue
        dp = [r[1] - r[3] for r in B]; dm = [r[2] - r[4] for r in B]
        sd = lambda x: (statistics.stdev(x) / math.sqrt(len(x))) if len(x) > 1 else float('nan')
        print(f'   [{lo},{hi}) n={len(B):3d}  p+1: {statistics.mean(dp):+6.1f} (se {sd(dp):.1f})   p-1: {statistics.mean(dm):+6.1f} (se {sd(dm):.1f})   bits {statistics.mean(r[0] for r in B):.0f}')
# fcy split by eps (D = 5 rows give eps = (5/p)); free Fermat side p - 1 is a different group from the forced p - eps only when inert
eps5 = {}
for line in open(os.path.join(d, 'fcy.rows')):
    t = line.split()
    if t[0] == 'fcy' and int(t[3]) == 5: eps5[int(t[1])] = int(t[4])
for e in (-1, 1):
    R = [(int(t[3]), float(t[4]), float(t[5]), float(t[6]), float(t[7])) for t in rows if t[1] == 'fcy' and eps5.get(int(t[2])) == e]
    print(f' fcy eps={e}: n={len(R)} (free side = p-1 only when eps = -1)')
    for lo, hi in ((0, 128), (128, 256), (256, 512), (512, 1024), (1024, 4096)):
        B = [r for r in R if lo <= r[0] < hi]
        if not B: continue
        dp = [r[1] - r[3] for r in B]; dm = [r[2] - r[4] for r in B]
        print(f'   [{lo},{hi}) n={len(B):3d}  p+1: {statistics.mean(dp):+6.1f}   p-1: {statistics.mean(dm):+6.1f}   mean bits {statistics.mean(r[0] for r in B):.0f}  (p-1 excess / bits {statistics.mean(dm)/statistics.mean(r[0] for r in B):.3f})')
nomatch = [t for t in rows if t[1] == 'fcy' and int(t[2]) not in eps5]
print(' fcy INH rows without a D=5 row:', len(nomatch))
