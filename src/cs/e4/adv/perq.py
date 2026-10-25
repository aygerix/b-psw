#!/usr/bin/env python3
"""per-q breakdown of small-q index hits (q | idx among records with q | No), from adv.py-verified inputs.
usage: perq.py FILE...   (PHI/REC/CTL files; recomputes orders like adv.py)"""
import sys
sys.argv = [sys.argv[0]] + sys.argv[1:]
import adv
phi, rec, ctl, _ = adv.load(sys.argv[1:])
fs, _ = adv.forced_set(phi, rec, verify=False)
Q = [3, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
from collections import defaultdict
cnt = defaultdict(lambda: [0, 0])
def add(tag, N, fa, o):
    No, oo = adv.oddpart(N), adv.oddpart(o)
    idx = No // oo
    for q in Q:
        if No % q == 0:
            cnt[(tag, q)][0] += 1
            cnt[(tag, q)][1] += (idx % q == 0)
for mode, x, p in fs:
    fa = rec[(mode, x, p)]
    if fa is None:
        continue
    if mode == 'F':
        N = p + 1; o = adv.order_from(N, fa, lambda m: adv.fibmod(m, p) == 0); add('FF', N, fa, o)
    else:
        N = p - 1; o = adv.order_from(N, fa, lambda m: pow(2, m, p) == 1); add('LF', N, fa, o)
for tag, x, p, fa in ctl:
    if fa is None:
        continue
    if tag == 'FFc':
        N = p + 1; o = adv.order_from(N, fa, lambda m: adv.fibmod(m, p) == 0)
    else:
        N = p - 1; o = adv.order_from(N, fa, lambda m: pow(2, m, p) == 1)
    add(tag, N, fa, o)
for tag in ('FF', 'FFc', 'LF', 'LFc'):
    tot_n = tot_h = 0; tot_e = 0.0
    parts = []
    for q in Q:
        n, h = cnt[(tag, q)]
        if n:
            parts.append('q%d %d/%d (exp %.1f)' % (q, h, n, n / q))
            tot_h += h; tot_e += n / q
    print(tag, ' '.join(parts[:6]), ' | total hits %d exp %.1f' % (tot_h, tot_e))
