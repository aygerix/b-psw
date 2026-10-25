#!/usr/bin/env python3
"""Independent summary of the verified E4 records (pickles from verify.py) + targeted checks.
Usage: summ.py   (expects verified_aF.pkl, verified_aL.pkl, verified_apilR_arep1_arep2.pkl in cwd)
"""
import math
import pickle
import random
from sympy import factorint, isprime, jacobi_symbol, primerange

TT = [3, 10, 100, 1000]
TU = [200, 1500, 10000]


def load(*names):
    out = []
    for n in names:
        out += pickle.load(open(n, 'rb'))
    return out


def z_(o, p, v):
    return (o - p) / math.sqrt(v) if v > 0 else float('nan')


def summ(rs):
    if not rs:
        return None
    d = {'n': len(rs)}
    so, sp, sv = sum(r['L'] for r in rs), sum(r['EL'] for r in rs), sum(r['VL'] for r in rs)
    d['log'] = (so, sp, z_(so, sp, sv))
    ho, he = sum(r['hq'] for r in rs), sum(r['ehq'] for r in rs)
    d['hq'] = (ho, he, z_(ho, he, he))
    for k, t in enumerate(TT):
        o = sum(r['o'][k] for r in rs); p = sum(r['p'][k] for r in rs)
        d['t%d' % t] = (o, p, z_(o, p, sum(r['p'][k] * (1 - r['p'][k]) for r in rs)))
    for k, T in enumerate(TU):
        o = sum(r['U'][k] for r in rs); p = sum(r['P'][k] for r in rs)
        d['U%d' % T] = (o, p, z_(o, p, sum(r['P'][k] * (1 - r['P'][k]) for r in rs)))
    return d


STRATA = [('all', lambda r: True), ('p<2^33', lambda r: r['bp'] <= 33), ('p>2^33', lambda r: r['bp'] > 33),
          ('p>2^64', lambda r: r['bp'] > 64)]
KEYS = ['log', 'hq', 't3', 't10', 't100', 't1000', 'U200', 'U1500', 'U10000']


def table(title, rows, side):
    ctl = side + 'c'
    print('=' * 90)
    print(title, side)
    for sn, sel in STRATA:
        a = summ([r for r in rows if r['tag'] == side and sel(r)])
        b = summ([r for r in rows if r['tag'] == ctl and sel(r)])
        if not a:
            continue
        print(' -- %s n=%d nc=%d' % (sn, a['n'], b['n'] if b else 0))
        for k in KEYS:
            x = a[k]
            s = '   %-6s obs %8.2f pred %8.2f z %+5.2f ratio %6.3f' % (k, x[0], x[1], x[2], x[0] / x[1] if x[1] else float('nan'))
            if b:
                y = b[k]
                rb = y[0] / y[1] if y[1] else float('nan')
                s += ' | ctl obs %8.2f pred %8.2f z %+5.2f ratio %6.3f' % (y[0], y[1], y[2], rb)
                ra = x[0] / x[1] if x[1] else float('nan')
                s += ' | kappa %6.3f' % (ra / rb if rb and rb == rb else float('nan'))
            print(s)


def v2(n):
    return (n & -n).bit_length() - 1


def gpf_odd(n):
    f = factorint(n)
    qs = [q for q in f if q != 2]
    return max(qs) if qs else 1


def fac_str(n):
    f = factorint(n)
    return '*'.join('%d^%d' % (q, a) if a > 1 else str(q) for q, a in sorted(f.items()))


def orders(r):
    if r['tag'] in ('FF', 'FFc'):
        return (r['ef'] if r['tag'] == 'FF' else None, r['ordr']) if r['tag'] == 'FF' else (None, r['ordr'])
    return (r['ordr'], r['ef'] if r['tag'] == 'LF' else None)


def main():
    main_ = load('verified_aF.pkl', 'verified_aL.pkl')
    rep = load('verified_apilR_arep1_arep2.pkl')
    rep_np = [r for r in rep if r['run'] != 'apilR']
    comb = main_ + rep
    table('MAIN', main_, 'FF')
    table('MAIN', main_, 'LF')
    table('REPLICATION (pilot + rep1 + rep2)', rep, 'LF')
    table('REPLICATION WITHOUT PILOT (rep1 + rep2)', rep_np, 'LF')
    table('COMBINED LF', comb, 'LF')

    # ---- usable forced primes p > 2^33: free side odd order T-smooth AND forced order odd part T-smooth
    print('=' * 90)
    print('Forced primes p > 2^33 usable on BOTH sides at T (free side odd order and forced order odd part T-smooth)')
    for T in TU:
        k = TU.index(T)
        lst = [r for r in comb if r['tag'] in ('FF', 'LF') and r['bp'] > 33 and r['U'][k] == 1 and gpf_odd(r['ef']) <= T]
        print(' T=%d: %d primes' % (T, len(lst)))
        for r in lst:
            p = r['p']
            if r['tag'] == 'FF':
                o2, z = r['ef'], r['ordr']
            else:
                o2, z = r['ordr'], r['ef']
            # independent recheck of both orders from scratch with sympy factorint
            assert pow(2, o2, p) == 1 and all(pow(2, o2 // q, p) != 1 for q in factorint(o2))
            print('   %s %s=%d p=%d bits=%d p%%8=%d ord2=%s z=%s class=(%d,%d) run=%s'
                  % (r['tag'], 'e' if r['tag'] == 'FF' else 'f', r['ef'], p, r['bp'], p % 8, fac_str(o2), fac_str(z),
                     v2(o2), v2(z), r['run']))
    # ---- main-run LF p > 2^33 with odd index >= 10
    print('=' * 90)
    print('MAIN LF p > 2^33 with odd Fermat index >= 10:')
    for r in main_:
        if r['tag'] == 'LF' and r['bp'] > 33 and r['idx'] >= 10:
            g = math.gcd(r['idx'], r['ef'])
            print('   f=%d (%s) bits=%d idx=%s gcd(idx,f)=%d idx mod f=%d' % (r['ef'], fac_str(r['ef']), r['bp'], fac_str(r['idx']), g, r['idx'] % r['ef']))
    print('REPLICATION LF p > 2^33 with odd Fermat index >= 10:')
    for r in rep:
        if r['tag'] == 'LF' and r['bp'] > 33 and r['idx'] >= 10:
            print('   f=%d bits=%d idx=%s run=%s' % (r['ef'], r['bp'], fac_str(r['idx']), r['run']))
    print('Max LF p bits: main %d, replication %d; max FF p bits %d' % (
        max(r['bp'] for r in main_ if r['tag'] == 'LF'), max(r['bp'] for r in rep if r['tag'] == 'LF'),
        max(r['bp'] for r in main_ if r['tag'] == 'FF')))
    print('Counts: FF %d (p>2^33 %d, p>2^64 %d), FFc %d; LF main %d rep %d (p>2^33 %d); LFc %d' % (
        sum(r['tag'] == 'FF' for r in main_), sum(r['tag'] == 'FF' and r['bp'] > 33 for r in main_),
        sum(r['tag'] == 'FF' and r['bp'] > 64 for r in main_), sum(r['tag'] == 'FFc' for r in main_),
        sum(r['tag'] == 'LF' for r in main_), sum(r['tag'] == 'LF' for r in rep),
        sum(r['tag'] == 'LF' and r['bp'] > 33 for r in comb), sum(r['tag'] == 'LFc' for r in comb)))
    # ---- forced primes without any control
    for side in ('FF', 'LF'):
        n_nc = 0
        for i, r in enumerate(comb):
            if r['tag'] == side:
                nxt = comb[i + 1] if i + 1 < len(comb) else None
                if not (nxt and nxt['tag'] == side + 'c' and nxt['p'] != r['p']):
                    n_nc += 1
        print('%s forced primes with no control: %d' % (side, n_nc))

    # ---- order channel: precise unconditional null for forced primes of 34..80 bits
    print('=' * 90)
    print('ORDER CHANNEL: E[P_T] over K random inert primes of the same bit length and class (trial division by q <= 1e4)')
    random.seed(20261025)
    small = list(primerange(3, 10001))
    K = 400
    for side, pool in (('FF', main_), ('LF', comb)):
        exp = [0.0, 0.0, 0.0]; obs = [0, 0, 0]; n = 0; lst = []
        for r in pool:
            if r['tag'] != side or not (33 < r['bp'] <= 80):
                continue
            n += 1
            B, m = r['bp'], r['ef']
            mod = m * 5
            res = [c for c in range(mod) if (c % m == (1 if side == 'FF' else m - 1)) and c % 5 in (2, 3)]
            acc = [0.0, 0.0, 0.0]; got = 0; tries = 0
            while got < K and tries < 200000:
                tries += 1
                x = random.randrange(1 << (B - 1), 1 << B)
                x = x - x % mod + random.choice(res)
                if x.bit_length() != B or not isprime(x):
                    continue
                N = x + 1 if side == 'FF' else x - 1
                No = N >> v2(N)
                rem = No; sT = [1, 1, 1]
                for q in small:
                    if rem % q:
                        continue
                    pe = 1
                    while rem % q == 0:
                        rem //= q; pe *= q
                    for j, T in enumerate(TU):
                        if q <= T:
                            sT[j] *= pe
                for j in range(3):
                    acc[j] += sT[j] / No
                got += 1
            for j in range(3):
                exp[j] += acc[j] / max(got, 1)
                ok = r['U'][j] == 1 and gpf_odd(r['ef']) <= TU[j]
                obs[j] += ok
            if r['U'][1] == 1:
                lst.append((r['p'], r['ef'], acc[1] / max(got, 1)))
        print(' %s n=%d (34..80 bits): ' % (side, n) + '  '.join(
            'T=%d obs %d exp %.3f' % (TU[j], obs[j], exp[j]) for j in range(3)))
        for p, ef, e1 in lst:
            print('    usable@1500: p=%d ef=%d  null P for its size/class %.4f' % (p, ef, e1))


if __name__ == '__main__':
    main()
