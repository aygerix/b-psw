#!/usr/bin/env python3
"""Free-side (Lucas) unabsorbed-rough slope for Fermat-forced families, split vs inert members, >= 64 bits, per D.
oc_summary.py reports inert members only; for split members the Lucas group is p - 1 (shared with the forced side)."""
import sys
for fam, src in (('cyc', 'cyc'), ('wag', 'wag'), ('gm+', 'gm'), ('gm-', 'gm')):
    for D in (5, -7, -11, 13, -15, 17, -19, -23):
        R = [l.split() for l in open('rows/' + src + '.rows')]
        R = [r for r in R if r[0] == fam and int(r[3]) == D and int(r[2]) >= 64]
        for e in (-1, 1):
            S = [r for r in R if int(r[4]) == e]
            if not S: continue
            sl = sum(float(r[12]) for r in S) / sum(int(r[2]) for r in S)
            print(f'{fam} D={D} eps={e:+d} n={len(S)} free(Lucas) slope {sl:.3f} usable {sum(1 for r in S if int(r[9]) > 0 and int(r[10]) > 0)}')
