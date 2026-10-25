#!/usr/bin/env python3
"""E4 forced-index coupling (FIC): summarise fic.gp records.

Record (23 fields): tag ef bits(p) bits(No) logidx Elog Vlog hq ehq o3 o10 o100 o1000 p3 p10 p100 p1000
                    U200 P200 U1500 P1500 U1e4 P1e4
For each tag (FF, FFc, LF, LFc) and size stratum, compare observed with the per-prime Kummer prediction
(uniform element of the cyclic group of odd order No) and report kappa = (obs/pred)_forced / (obs/pred)_control.
Usage: analyze.py rec1 [rec2 ...]
"""
import math
import sys

TAILS = [3, 10, 100, 1000]
TU = [200, 1500, 10000]


def load(paths):
    rows = []
    for path in paths:
        for line in open(path):
            if line.startswith('#') or not line.strip():
                continue
            f = line.split()
            if len(f) != 23:
                continue
            r = dict(tag=f[0], ef=int(f[1]), bp=int(f[2]), bN=int(f[3]),
                     L=float(f[4]), EL=float(f[5]), VL=float(f[6]), hq=int(f[7]), ehq=float(f[8]))
            r['o'] = [int(x) for x in f[9:13]]
            r['p'] = [float(x) for x in f[13:17]]
            r['U'] = [int(f[17]), int(f[19]), int(f[21])]
            r['P'] = [float(f[18]), float(f[20]), float(f[22])]
            rows.append(r)
    return rows


def summ(rs):
    n = len(rs)
    if n == 0:
        return None
    out = {'n': n}
    so, sp, sv = sum(r['L'] for r in rs), sum(r['EL'] for r in rs), sum(r['VL'] for r in rs)
    out['log'] = (so, sp, (so - sp) / math.sqrt(sv) if sv > 0 else float('nan'))
    ho, he = sum(r['hq'] for r in rs), sum(r['ehq'] for r in rs)
    out['hq'] = (ho, he, (ho - he) / math.sqrt(he) if he > 0 else float('nan'))
    for k, t in enumerate(TAILS):
        o = sum(r['o'][k] for r in rs)
        p = sum(r['p'][k] for r in rs)
        v = sum(r['p'][k] * (1 - r['p'][k]) for r in rs)
        out['t%d' % t] = (o, p, (o - p) / math.sqrt(v) if v > 0 else float('nan'))
    for k, T in enumerate(TU):
        o = sum(r['U'][k] for r in rs)
        p = sum(r['P'][k] for r in rs)
        v = sum(r['P'][k] * (1 - r['P'][k]) for r in rs)
        out['U%d' % T] = (o, p, (o - p) / math.sqrt(v) if v > 0 else float('nan'))
    return out


def ratio(x):
    o, p = x[0], x[1]
    return o / p if p > 0 else float('nan')


def main():
    rows = load(sys.argv[1:])
    strata = [('all', lambda r: True), ('p<2^33', lambda r: r['bp'] <= 33),
              ('p>2^33', lambda r: r['bp'] > 33), ('p>2^64', lambda r: r['bp'] > 64)]
    keys = ['log', 'hq'] + ['t%d' % t for t in TAILS] + ['U%d' % T for T in TU]
    for side, ctl in (('FF', 'FFc'), ('LF', 'LFc')):
        print('=' * 100)
        print('%s (forced) vs %s (matched control): obs, pred, z, obs/pred; kappa = ratio_forced / ratio_control'
              % (side, ctl))
        for sname, sel in strata:
            a = summ([r for r in rows if r['tag'] == side and sel(r)])
            b = summ([r for r in rows if r['tag'] == ctl and sel(r)])
            if a is None:
                continue
            print('-- stratum %s: n_forced = %d, n_control = %d' % (sname, a['n'], b['n'] if b else 0))
            for k in keys:
                fa = a[k]
                line = '   %-7s forced obs %10.2f pred %10.2f z %+6.2f ratio %6.3f' % (k, fa[0], fa[1], fa[2], ratio(fa))
                if b:
                    fb = b[k]
                    line += ' | control obs %10.2f pred %10.2f z %+6.2f ratio %6.3f' % (fb[0], fb[1], fb[2], ratio(fb))
                    rb = ratio(fb)
                    line += ' | kappa %6.3f' % (ratio(fa) / rb if rb and rb == rb and rb > 0 else float('nan'))
                print(line)
    # source yield outside the census range (p > 1e10 ~ 2^33.2): free primes added to a pool
    print('=' * 100)
    print('Source yield for p > 2^34 (outside census2 range), forced order e (or f) with odd part T-smooth:')
    print('  free usable primes at rows <= T: count of U_T (obs) vs sum of P_T (Kummer pred)')

    def gpf_odd(m):
        while m % 2 == 0:
            m //= 2
        g, d = 1, 3
        while d * d <= m:
            while m % d == 0:
                g, m = d, m // d
            d += 2
        return max(g, m) if m > 1 else g

    for side in ('FF', 'LF'):
        rs = [r for r in rows if r['tag'] == side and r['bp'] >= 34]
        parts = []
        for k, T in enumerate(TU):
            sel = [r for r in rs if gpf_odd(r['ef']) <= T]
            parts.append('T=%d (n=%d) obs %d pred %.4g' % (T, len(sel), sum(r['U'][k] for r in sel),
                                                           sum(r['P'][k] for r in sel)))
        print('  %s n=%d  ' % (side, len(rs)) + '  '.join(parts))
    # unconditional yield: forced usability rate vs the matched control's PREDICTED rate (the control's group
    # orders are random-prime group orders, so this also tests coupling through the group order p -+ 1 itself)
    print('Unconditional yield per prime, p > 2^33: forced obs rate vs control pred rate (kappa_total)')
    for side, ctl in (('FF', 'FFc'), ('LF', 'LFc')):
        a = [r for r in rows if r['tag'] == side and r['bp'] > 33]
        b = [r for r in rows if r['tag'] == ctl and r['bp'] > 33]
        if not a or not b:
            continue
        parts = []
        for k, T in enumerate(TU):
            ra = sum(r['U'][k] for r in a) / len(a)
            rb = sum(r['P'][k] for r in b) / len(b)
            ra2 = sum(r['P'][k] for r in a) / len(a)
            parts.append('T=%d forced %.4g (pred-from-own-order %.4g) control %.4g kappa_total %.3g'
                         % (T, ra, ra2, rb, ra / rb if rb > 0 else float('nan')))
        print('  %s: n=%d, control n=%d\n    ' % (side, len(a), len(b)) + '\n    '.join(parts))
    bp = {}
    for r in rows:
        if r['tag'] in ('FF', 'LF'):
            bp.setdefault(r['tag'], []).append(r['bp'])
    for tag, v in bp.items():
        v.sort()
        print('  %s size (bits of p): min %d median %d max %d' % (tag, v[0], v[len(v) // 2], v[-1]))


if __name__ == '__main__':
    main()
