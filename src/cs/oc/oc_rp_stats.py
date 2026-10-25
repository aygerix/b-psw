#!/usr/bin/env python3
"""Riesel c 2^m - 1 / Proth c 2^m + 1 (odd c <= 127, m <= 600): per (D, class) member counts and usable counts.
Free side: Riesel -> Lucas (p+1 = c 2^m) when inert; Proth -> Fermat (p-1 = c 2^m). Other side generic."""
import sys
from collections import defaultdict
sys.path.insert(0, '.')
from oc_post import load
for fam in ('rie', 'pro'):
    R = load([f'out/full_{fam}.out'])
    print(f'=== {fam}: {len(set((r["fam"], r["par"]) for r in R))} members')
    for D in (5, -7, -11, 13, -15):
        RD = [r for r in R if r['D'] == D]
        cls = defaultdict(list)
        for r in RD: cls[(r['j1'], r['j2'])].append(r)
        rows = []
        for c, v in cls.items():
            u = [r for r in v if r['oo2'] > 0 and r['oz2'] > 0]
            rows.append((len(u), len(v), c, max([r['bits'] for r in u], default=0)))
        rows.sort(reverse=True)
        print(f'  D={D}: classes={len(cls)}; top by usable: ' + '; '.join(f'{c}: {nu}/{nv} (max usable bits {mb})' for nu, nv, c, mb in rows[:6]))
        # free-side boundedness check and other-side rough bits by size band
        free = 'rL' if fam == 'rie' else 'rF'; other = 'rF' if fam == 'rie' else 'rL'
        inert = [r for r in RD if r['eps'] == -1]
        for lo, hi in ((0, 40), (40, 100), (100, 200), (200, 400), (400, 800)):
            band = [r for r in inert if lo <= r['bits'] < hi]
            if not band: continue
            nf0 = sum(1 for r in band if r[free] == 0)
            avg_other = sum(r[other] for r in band) / len(band)
            nus = sum(1 for r in band if r['oo2'] > 0 and r['oz2'] > 0)
            print(f'     inert bits [{lo},{hi}): n={len(band)} free-side unabsorbed=0 for {nf0}; other side mean rough bits {avg_other:.1f} ({avg_other/((lo+hi)/2):.2f}/bit); usable {nus}')
