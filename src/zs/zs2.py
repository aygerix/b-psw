#!/usr/bin/env python3
"""DFS search for Psi-families (data from vec2.gp).

Korselt mode:  lcm(Psi_d - 1) | n - 1                      (Carmichael identity, row C3F)
Williams mode: also lcm(Psi_d + 1) | n + 1                  (rows C3F + C3L)
At each special point theta_e (root of some Psi_d -+ 1 in S): prod_S Psi_d(theta_e) = +1 (minus side) or -1
(plus side), i.e. unit exponents sum to 0 and torsion bit sums to 0 (resp. 1); where some member has
multiplicity 2, also sum_S Psi_d'/Psi_d(theta_e) = 0. Multiplicity >= 3 is reported, not tested.
Every reported family is re-verified as an exact polynomial identity in GP (verify.gp).
Usage: python3 zs2.py data.out korselt|williams [--kmax K] [--no5]
"""
import sys
from fractions import Fraction

def parse(path):
    C = []
    for line in open(path):
        f = line.strip().split('|')
        if len(f) < 5 or not f[0].isdigit():
            continue
        def sp(s):
            if s == 'X':
                return None
            out = {}
            for t in s.split(','):
                e, m = t.split(':'); out[int(e)] = int(m)
            return out
        st = {}
        for t in f[4:]:
            e, v = t.split('=', 1); e = int(e)
            if v in ('N', 'Z'):
                st[e] = None; continue
            u, dv = v[1:].split(';')
            if u.startswith('['):
                ex = [int(x) for x in u.strip('[]').split(',')]
                unit = tuple(ex[:-1]) + (ex[-1] % 2,)
                der = tuple(Fraction(x) for x in dv.strip('[]').split(','))
            else:
                unit = ((1 if int(u) < 0 else 0),)
                der = (Fraction(dv),)
            st[e] = (unit, der)
        C.append(dict(d=int(f[0]), deg=int(f[1]), m=sp(f[2]), p=sp(f[3]), st=st))
    return C

def add(a, b):
    return tuple(x + y for x, y in zip(a, b))

def check(S, C, Tm, Tp):
    """exact zero-sum conditions for the set S (indices) with special sets Tm, Tp (dicts e -> max mult)."""
    for side, T in (('m', Tm), ('p', Tp)):
        for e, mult in T.items():
            if mult >= 3:
                return 'mult3'
            us = None; ds = None
            for i in S:
                u, dv = C[i]['st'][e]
                us = u if us is None else add(us, u)
                ds = dv if ds is None else add(ds, dv)
            if any(us[:-1]) or us[-1] % 2 != (0 if side == 'm' else 1):
                return False
            if mult == 2 and any(ds):
                return False
    return True

def search(C, mode, kmax, no5):
    W = mode == 'williams'
    idx = [i for i, c in enumerate(C) if c['m'] is not None and (not W or c['p'] is not None)
           and (not no5 or c['d'] % 5)]
    found, mult3, nodes = [], [], [0]

    def rec(start, S, Tm, Tp):
        nodes[0] += 1
        if len(S) >= 2:
            r = check(S, C, Tm, Tp)
            if r is True:
                found.append(tuple(C[i]['d'] for i in S))
            elif r == 'mult3':
                mult3.append(tuple(C[i]['d'] for i in S))
        if len(S) == kmax:
            return
        for j in idx:
            if j <= start:
                continue
            c = C[j]
            Tm2 = dict(Tm)
            for e, mm in c['m'].items():
                Tm2[e] = max(Tm2.get(e, 0), mm)
            Tp2 = dict(Tp)
            if W:
                for e, mm in c['p'].items():
                    Tp2[e] = max(Tp2.get(e, 0), mm)
                if set(Tm2) & set(Tp2):
                    continue
            allT = set(Tm2) | set(Tp2)
            ok = all(C[i]['st'].get(e) is not None for i in S + [j] for e in allT)
            if ok:
                rec(j, S + [j], Tm2, Tp2)

    rec(-1, [], {}, {})
    return found, mult3, nodes[0]

if __name__ == '__main__':
    C = parse(sys.argv[1]); mode = sys.argv[2]
    kmax = int(sys.argv[sys.argv.index('--kmax') + 1]) if '--kmax' in sys.argv else 8
    no5 = '--no5' in sys.argv
    F, M3, nodes = search(C, mode, kmax, no5)
    print(f"mode {mode}, candidates {len(C)}, kmax {kmax}{', 5 not | d' if no5 else ''}: "
          f"DFS nodes {nodes}, families {len(F)}, untested (multiplicity >= 3) {len(M3)}")
    for s in sorted(F, key=lambda s: (len(s), s)):
        print('FAMILY', ' '.join(map(str, s)))
    for s in M3[:20]:
        print('MULT3', ' '.join(map(str, s)))
