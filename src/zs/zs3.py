#!/usr/bin/env python3
"""Mixed Psi-families: every member satisfies the Korselt side (Psi_d - 1 | n - 1 via special points);
members of type W also satisfy Psi_d + 1 | n + 1 (Williams side, n(theta) = -1 at their plus-special points);
members of type F are Lucas-forced at Y = L_2K (D = 5) and need 5 not | d (Lemma 12(b)).
Requires >= 1 member of each type. Uses parse/check helpers of zs2.py.  Usage: python3 zs3.py data.out [--kmax K]"""
import sys
from zs2 import parse, add

def check(S, C, Tm, Tp):
    for side, T in (('m', Tm), ('p', Tp)):
        for e, mult in T.items():
            if mult >= 3:
                return 'mult3'
            us = ds = None
            for i, _ in S:
                u, dv = C[i]['st'][e]
                us = u if us is None else add(us, u); ds = dv if ds is None else add(ds, dv)
            if any(us[:-1]) or us[-1] % 2 != (0 if side == 'm' else 1):
                return False
            if mult == 2 and any(ds):
                return False
    return True

def main():
    C = parse(sys.argv[1])
    kmax = int(sys.argv[sys.argv.index('--kmax') + 1]) if '--kmax' in sys.argv else 8
    opts = []
    for i, c in enumerate(C):
        if c['m'] is None:
            continue
        if c['d'] % 5:
            opts.append((i, 'F'))
        if c['p'] is not None:
            opts.append((i, 'W'))
    found, nodes = [], [0]
    def rec(start, S, Tm, Tp):
        nodes[0] += 1
        types = {t for _, t in S}
        if len(S) >= 2 and types == {'F', 'W'} and check(S, C, Tm, Tp) is True:
            found.append(tuple((C[i]['d'], t) for i, t in S))
        if len(S) == kmax:
            return
        used = {i for i, _ in S}
        for j, (i, t) in enumerate(opts):
            if j <= start or i in used:
                continue
            Tm2 = dict(Tm)
            for e, mm in C[i]['m'].items():
                Tm2[e] = max(Tm2.get(e, 0), mm)
            Tp2 = dict(Tp)
            if t == 'W':
                for e, mm in C[i]['p'].items():
                    Tp2[e] = max(Tp2.get(e, 0), mm)
            if set(Tm2) & set(Tp2):
                continue
            allT = set(Tm2) | set(Tp2)
            if all(C[a]['st'].get(e) is not None for a in [x for x, _ in S] + [i] for e in allT):
                rec(j, S + [(i, t)], Tm2, Tp2)
    rec(-1, [], {}, {})
    print(f"mixed search, candidates {len(C)}, options {len(opts)}, kmax {kmax}: DFS nodes {nodes[0]}, families {len(found)}")
    for f in found:
        print('FAMILY', ' '.join(f"{d}{t}" for d, t in f))

if __name__ == '__main__':
    main()
