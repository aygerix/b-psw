#!/usr/bin/env python3
"""Line OC post-processing.
Input rows (oc_fam.gp): fam par bits D eps j1 j2 oo200 oz200 oo1500 oz1500 roughF roughL
  oo/oz = exact odd(ord_p 2) / odd(z(p)) when B-smooth, else 0.
For every (D, class (j1,j2)) list members, usable ones (both odd orders 1500-smooth), and run the exact
subset search on usable members: Lemma 1 criterion (one class, odd orders | n -+ 1, Selfridge picks D,
v2 room), then src/bpsw.py on survivors. Members are regenerated from (fam, par) exactly.
"""
import sys, itertools, math, os
from collections import defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from bpsw import bpsw_pass, jacobi, selfridge_D

def fib(n):
    a, b = 0, 1
    for _ in range(n): a, b = b, a + b
    return a
def member(fam, par):
    if fam == 'mer': return 2**par - 1
    if fam == 'wag': return (2**par + 1) // 3
    if fam == 'fib': return fib(par)
    if fam == 'luc': return fib(par - 1) + fib(par + 1)
    if fam == 'gm+': k = (par - 1)//2; return 2**par + 2**(k+1) + 1
    if fam == 'gm-': k = (par - 1)//2; return 2**par - 2**(k+1) + 1
    if fam.startswith('rie'): return int(fam[3:]) * 2**par - 1
    if fam.startswith('pro'): return int(fam[3:]) * 2**par + 1
    if fam == 'cyc':
        from sympy import cyclotomic_poly
        from math import gcd
        P = int(cyclotomic_poly(par, 2)); return P // gcd(P, par)
    if fam == 'fcy':
        from sympy import divisors
        from sympy.ntheory import mobius
        from math import gcd
        num = den = 1
        for d in divisors(par):
            mu = mobius(par // d)
            if mu == 1: num *= fib(d)
            elif mu == -1: den *= fib(d)
        assert num % den == 0
        P = num // den; return P // gcd(P, par)
    raise ValueError(fam)

FORCED = {'fcy': 'L', 'mer': 'F', 'wag': 'F', 'gm+': 'F', 'gm-': 'F', 'cyc': 'F', 'fib': 'L', 'luc': 'L'}
def forced_side(fam):
    if fam.startswith('rie'): return 'L*'   # Lucas order bounded pieces (p+1 = c 2^m), Fermat generic
    if fam.startswith('pro'): return 'F*'   # Fermat order bounded pieces (p-1 = c 2^m), Lucas generic
    return FORCED[fam]

def load(paths):
    rows = []
    for path in paths:
        for line in open(path):
            t = line.split()
            if len(t) != 13: continue
            r = dict(fam=t[0], par=int(t[1]), bits=int(t[2]), D=int(t[3]), eps=int(t[4]), j1=int(t[5]), j2=int(t[6]),
                     oo2=int(t[9]), oz2=int(t[10]), rF=float(t[11]), rL=float(t[12]))
            rows.append(r)
    return rows

def selfridge_seq(Dstar):
    """All Selfridge D values strictly before Dstar."""
    out, D = [], 5
    while D != Dstar:
        out.append(D); D = -D - 2 if D > 0 else -D + 2
        if abs(D) > 1000: raise ValueError
    return out

def subset_search(mem, D, j1, j2, maxk=None, log=print):
    """mem: list of (p, oo, oz, eps). Exact Lemma-1 prefilter then bpsw.py."""
    before = selfridge_seq(D)
    hits, tested, passed_pref = [], 0, 0
    m = len(mem)
    maxk = m if maxk is None else min(m, maxk)
    for k in range(2, maxk + 1):
        for S in itertools.combinations(range(m), k):
            tested += 1
            # Jacobi: (D/n) = prod eps = -1
            e = 1
            for i in S: e *= mem[i][3]
            if e != -1: continue
            L1 = 1; L2 = 1
            for i in S:
                L1 = L1 * mem[i][1] // math.gcd(L1, mem[i][1]); L2 = L2 * mem[i][2] // math.gcd(L2, mem[i][2])
            if math.gcd(L1, L2) != 1: continue
            n = 1
            for i in S: n *= mem[i][0]
            if (n - 1) % L1 or (n + 1) % L2: continue
            if any(jacobi(Dp % n if Dp > 0 else Dp % n, n) == -1 for Dp in before): continue
            s = ((n - 1) & -(n - 1)).bit_length() - 1; t = ((n + 1) & -(n + 1)).bit_length() - 1
            if j1 > s or j2 > t: continue
            passed_pref += 1
            ok = bpsw_pass(n)
            log(f'  PREFILTER PASS D={D} class=({j1},{j2}) S={[mem[i][0] if mem[i][0] < 10**12 else "~2^%d" % mem[i][0].bit_length() for i in S]} bpsw={ok} selD={selfridge_D(n)}')
            if ok: hits.append([mem[i][0] for i in S])
    return tested, passed_pref, hits

def main(paths, maxk=None):
    rows = load(paths)
    fams = sorted(set(r['fam'] if not (r['fam'].startswith('rie') or r['fam'].startswith('pro')) else r['fam'][:3] for r in rows))
    for F in fams:
        R = [r for r in rows if (r['fam'] == F or (F in ('rie', 'pro') and r['fam'].startswith(F)))]
        nmem = len(set((r['fam'], r['par']) for r in R))
        print(f'=== family {F}: {nmem} members, bits {min(r["bits"] for r in R)}..{max(r["bits"] for r in R)}')
        for D in sorted(set(r['D'] for r in R), key=lambda d: (abs(d), d)):
            RD = [r for r in R if r['D'] == D]
            ninert = sum(1 for r in RD if r['eps'] == -1)
            cls = defaultdict(list)
            for r in RD: cls[(r['j1'], r['j2'])].append(r)
            usable = {c: [r for r in v if r['oo2'] > 0 and r['oz2'] > 0] for c, v in cls.items()}
            big = sorted(cls.items(), key=lambda kv: -len(kv[1]))[:4]
            desc = ' '.join(f'({c[0]},{c[1]}):{len(v)}/u{len(usable[c])}' for c, v in big)
            print(f'  D={D:4d} members={len(RD)} inert={ninert} classes={len(cls)} top: {desc}')
            for c, U in usable.items():
                if len(U) >= 2:
                    mem = []
                    for r in U:
                        mem.append((member(r['fam'], r['par']), r['oo2'], r['oz2'], r['eps']))
                    if len(mem) >= 2:
                        tested, pp, hits = subset_search(mem, D, c[0], c[1], maxk)
                        print(f'    subset search D={D} class={c}: usable={len(mem)} subsets={tested} prefilter_pass={pp} BPSW_hits={len(hits)}')
                        for h in hits: print('    *** HIT', h)

if __name__ == '__main__':
    mk = int(os.environ.get('MAXK', '0')) or None
    main(sys.argv[1:], mk)
