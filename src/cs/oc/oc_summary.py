#!/usr/bin/env python3
"""Line OC summary table (D = 5 and D = -7). For each family:
 forced side, free (opposite) side, measured unabsorbed free-side rough bits at B = 1500 vs member bits (slope),
 achievable classes and the largest class, largest class restricted to members whose BOTH odd orders are 1500-smooth,
 per-member group-bit cost = log2(odd order on the forced side) + unabsorbed rough bits on the free side.
 Frontier reference: 0.16-0.22 members per group bit at small b, i.e. 4.5-6.3 bits per member."""
import sys, math, statistics
from collections import defaultdict
sys.path.insert(0, '.')
from oc_post import load
fams = [('mer', 'F', 'L'), ('wag', 'F', 'L'), ('gm+', 'F', 'L'), ('gm-', 'F', 'L'), ('cyc', 'F', 'L'),
        ('fib', 'L', 'F'), ('luc', 'L', 'F'), ('fcy', 'L', 'F')]
files = {'mer': 'out/full_mer.out', 'wag': 'out/full_wag.out', 'gm+': 'out/full_gm.out', 'gm-': 'out/full_gm.out',
         'cyc': 'out/full_cyc.out', 'fib': 'out/full_fib.out', 'luc': 'out/full_luc.out', 'fcy': 'out/full_fcy.out'}
def forced_order_bits(fam, par):
    if fam in ('mer', 'wag', 'gm+', 'gm-'): return math.log2(par)          # odd part of ord_p 2 = q
    if fam == 'cyc': e = par; o = e >> ((e & -e).bit_length() - 1); return math.log2(o) if o > 1 else 0.0
    if fam in ('fib', 'luc'): return math.log2(par)
    if fam == 'fcy': f = par; o = f >> ((f & -f).bit_length() - 1); return math.log2(o) if o > 1 else 0.0                          # odd part of z = f
for D in (5, -7):
    print(f'================ D = {D} ================')
    print(f'{"fam":4s} {"n":>3s} {"bits":>10s} {"inert":>5s} {"#cls":>4s} {"maxcls":>6s} {"maxcls_u":>8s} {"slope":>6s} {"cost_med(>=64b)":>15s} {"cost_min(>=64b)":>15s}')
    for fam, fs, fr in fams:
        R = [r for r in load([files[fam]]) if r['fam'] == fam and r['D'] == D]
        if not R: continue
        cls = defaultdict(list)
        for r in R: cls[(r['j1'], r['j2'])].append(r)
        maxc = max(len(v) for v in cls.values())
        maxu = max(sum(1 for r in v if r['oo2'] > 0 and r['oz2'] > 0) for v in cls.values())
        big = [r for r in R if r['bits'] >= 64]
        free = 'rL' if fr == 'L' else 'rF'
        # for split members the "free" Lucas group is p - 1, not the algebraic side: report inert-only slope for F-forced
        sel = [r for r in big if (fr == 'F' or r['eps'] == -1)]
        slope = (sum(r[free] for r in sel) / sum(r['bits'] for r in sel)) if sel else float('nan')
        # Lucas forcing (fib, luc, fcy) holds only for D = 5: for other D both sides are generic
        lforced_lost = (fs == 'L' and D != 5)
        costs = [(forced_order_bits(fam, r['par']) + r[free] + (r['rL'] if lforced_lost else 0.0)) if not (r['oo2'] > 0 and r['oz2'] > 0)
                 else math.log2(r['oo2']) + math.log2(r['oz2']) for r in sel]
        cm = statistics.median(costs) if costs else float('nan'); cmin = min(costs) if costs else float('nan')
        # exact odd-order bits for usable members (both odd orders 1500-smooth): upper bound on new group bits
        ub = [math.log2(r['oo2']) + math.log2(r['oz2']) for r in sel if r['oo2'] > 0 and r['oz2'] > 0]
        cmin = min(cmin, min(ub)) if ub and min(ub) > cmin else cmin
        print(f'{fam:4s} {len(R):3d} {min(r["bits"] for r in R):4d}-{max(r["bits"] for r in R):<5d} {sum(1 for r in R if r["eps"] == -1):5d} {len(cls):4d} {maxc:6d} {maxu:8d} {slope:6.2f} {cm:15.1f} {cmin:15.1f}')
    # Mersenne detail: inert members' j2 vs q
    if D in (5, -7):
        M = [r for r in load([files['mer']]) if r['D'] == D and r['eps'] == -1]
        print('  Mersenne inert (q, j2, q-1-j2):', [(r['par'], r['j2'], r['par'] - 1 - r['j2']) for r in M])
