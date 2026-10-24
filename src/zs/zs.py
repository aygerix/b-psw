#!/usr/bin/env python3
"""Zero-sum search for Korselt-identity families of Psi_d (src/zs/vec.out).

A set S of Psi_d satisfies lcm(Psi_d - 1) | prod Psi_d - 1 in Q[Y] iff, at every special point theta_e
(root of some Psi_d - 1, d in S): prod_S Psi_d(theta_e) = 1 (unit exponents and torsion sum to 0), and
sum_S Psi_d'/Psi_d(theta_e) = 0 where some member has multiplicity 2 (multiplicity >= 3 not handled).
Enumerates the exact special set T = E(S) and multiplicity set T2, then meet-in-the-middle per (T, T2).
Usage: python3 zs.py vec.out [--no5] [--kmax K]
"""
import sys, itertools
from fractions import Fraction
from collections import defaultdict

def parse(path):
    C = []
    for line in open(path):
        f = line.strip().split('|')
        if len(f) < 4:
            continue
        d, deg = int(f[0]), int(f[1])
        sp = {}
        if f[2]:
            for t in f[2].split(','):
                e, m = t.split(':'); sp[int(e)] = int(m)
        st = {}
        for t in f[3:]:
            e, v = t.split('=', 1); e = int(e)
            if v in ('N', 'Z'):
                st[e] = None; continue
            u, dv = v[1:].split(';')
            if u.startswith('['):
                ex = [int(x) for x in u.strip('[]').split(',')] if u.strip('[]') else []
                unit = tuple(ex[:-1]) + (ex[-1] % 2,)
                der = tuple(Fraction(x) for x in dv.strip('[]').split(','))
            else:
                unit = ((1 if int(u) < 0 else 0),)
                der = (Fraction(dv),)
            st[e] = (unit, der)
        C.append((d, deg, sp, st))
    return C

def vec(c, T, T2):
    d, deg, sp, st = c
    v = []
    for e in T:
        v += list(st[e][0])
    for e in T2:
        v += list(st[e][1])
    return v

def search(C, kmax, no5):
    E_all = sorted({e for c in C for e in c[3]})
    found = set()
    # candidate special sets: unions of members' special sets are what matter; enumerate T from members
    Ts = set()
    for c in C:
        Ts.add(frozenset(c[2]))
    # all unions of up to 3 members' special sets (families are unions of their members' sets)
    base = list(Ts)
    unions = set(base)
    for a in base:
        for b in base:
            unions.add(a | b)
    for a in list(unions):
        for b in base:
            unions.add(a | b)
    for T in sorted(unions, key=lambda s: (len(s), sorted(s))):
        T = sorted(T)
        cand = [c for c in C if set(c[2]) <= set(T) and all(c[3].get(e) is not None for e in T)
                and (not no5 or c[0] % 5) and max(c[2].values()) <= 2]
        if len(cand) < 2:
            continue
        mult_es = sorted({e for c in cand for e, m in c[2].items() if m >= 2})
        for r in range(len(mult_es) + 1):
            for T2 in itertools.combinations(mult_es, r):
                cc = [c for c in cand if all(m == 1 or e in T2 for e, m in c[2].items())]
                if len(cc) < 2:
                    continue
                V = [vec(c, T, T2) for c in cc]
                nt = sum(len(cc[0][3][e][0]) for e in T)  # unit coords; torsion bits sit at known places
                tors_idx = []
                pos = 0
                for e in T:
                    L = len(cc[0][3][e][0]); tors_idx.append(pos + L - 1); pos += L
                def key(s):
                    return tuple((x % 2) if i in tors_idx else x for i, x in enumerate(s))
                n = len(cc); h = n // 2
                A, B = list(range(h)), list(range(h, n))
                tab = defaultdict(list)
                for k in range(0, min(kmax, len(A)) + 1):
                    for sub in itertools.combinations(A, k):
                        s = [0] * len(V[0])
                        for i in sub:
                            s = [x + y for x, y in zip(s, V[i])]
                        tab[key(s)].append(sub)
                for k in range(0, min(kmax, len(B)) + 1):
                    for sub in itertools.combinations(B, k):
                        s = [0] * len(V[0])
                        for i in sub:
                            s = [-(x + y) if False else x + y for x, y in zip(s, V[i])]
                        neg = key([-x for x in s])
                        for sa in tab.get(neg, []):
                            S = sa + sub
                            if len(S) < 2 or len(S) > kmax:
                                continue
                            ds = tuple(sorted(cc[i][0] for i in S))
                            ES = set().union(*[set(cc[i][2]) for i in S])
                            M2 = {e for i in S for e, m in cc[i][2].items() if m >= 2}
                            if ES == set(T) and M2 == set(T2):
                                found.add(ds)
    return sorted(found, key=lambda s: (len(s), s))

if __name__ == '__main__':
    C = parse(sys.argv[1])
    no5 = '--no5' in sys.argv
    kmax = int(sys.argv[sys.argv.index('--kmax') + 1]) if '--kmax' in sys.argv else 8
    F = search(C, kmax, no5)
    print(f"candidates {len(C)}, families found (k <= {kmax}{', 5 not | d' if no5 else ''}): {len(F)}")
    for s in F:
        print(' '.join(map(str, s)))
