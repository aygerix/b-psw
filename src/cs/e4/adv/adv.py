#!/usr/bin/env python3
"""Adversarial re-check of E4 (FIC). Written independently of fic.gp / analyze.py.

  adv.py gen K SEED OUT.gp FILE...   draw K matched controls per forced prime (same bit length, p' == 1 mod e
                                     (FF) or p' == -1 mod f (LF), p' == 2,3 mod 5), write a gp vector CT.
  adv.py stats [--set NAME=FILE,FILE ...] FILE...  verify every record and print the statistics.

Every factorisation is checked (product + sympy.isprime on every factor). Phi_e(2) is recomputed with
sympy.cyclotomic_poly, Phi*_f as prod_{d|f} F_d^mu(f/d) (the gp side used the Moebius product for Phi_e(2) and the
homogenised cyclotomic value for Phi*_f, so the two sides use different formulas). Forcing (ord_p 2 = e exactly,
z(p) = f exactly), inertness and primitivity are re-derived here; the other order is computed here with a
Fibonacci fast-doubling / pow based order routine.

Kummer model (H0): for odd q with q^a || No (odd part of the group order), P(v_q(idx) >= j) = q^-j (0 <= j <= a),
independent over q. E[v] = sum_{j>=1} P(v>=j), E[v^2] = sum_{j>=1} (2j-1) P(v>=j).
"""
import math
import random
import sys
from collections import defaultdict

from sympy import isprime, jacobi_symbol, cyclotomic_poly, divisors, primefactors
from sympy.functions.combinatorial.numbers import mobius

TT = [3, 10, 100, 1000]
TU = [200, 1500, 10000]


def parse_fac(s):
    out = {}
    if s == '1':
        return out
    for t in s.split('*'):
        q, a = t.split('^')
        out[int(q)] = out.get(int(q), 0) + int(a)
    return out


def prod_fac(fa):
    r = 1
    for q, a in fa.items():
        r *= q ** a
    return r


def check_fac(N, fa, what):
    if prod_fac(fa) != N:
        raise SystemExit('BAD product %s' % what)
    for q in fa:
        if not isprime(q):
            raise SystemExit('BAD non-prime factor %d in %s' % (q, what))


def fib_pair(n, m):
    """(F_n mod m, F_{n+1} mod m) by fast doubling."""
    a, b = 0, 1
    for bit in bin(n)[2:]:
        c = (a * ((2 * b - a) % m)) % m
        d = (a * a + b * b) % m
        if bit == '1':
            a, b = d, (c + d) % m
        else:
            a, b = c, d
    return a, b


def fibmod(n, m):
    return fib_pair(n, m)[0]


def fib_exact(n):
    a, b = 0, 1
    for bit in bin(n)[2:]:
        c = a * (2 * b - a)
        d = a * a + b * b
        if bit == '1':
            a, b = d, c + d
        else:
            a, b = c, d
    return a


_fibcache = {}


def phistar(f):
    if f in _fibcache:
        return _fibcache[f]
    num, den = 1, 1
    for d in divisors(f):
        mu = int(mobius(f // d))
        if mu == 1:
            num *= fib_exact(d)
        elif mu == -1:
            den *= fib_exact(d)
    assert num % den == 0
    _fibcache[f] = num // den
    return num // den


def phi2(e):
    return int(cyclotomic_poly(e, 2))


def order_from(N, fa, is_one):
    """smallest divisor o of N with is_one(o), given N's factorisation (is_one(N) must hold)."""
    o = N
    for q, a in fa.items():
        for _ in range(a):
            if o % q == 0 and is_one(o // q):
                o //= q
            else:
                break
    return o


def v2(n):
    return (n & -n).bit_length() - 1


def oddpart(n):
    return n >> v2(n)


def kummer_stats(N, fa, o):
    """per-prime record: observed and predicted statistics for an element of order o in cyclic group of order N."""
    No, oo = oddpart(N), oddpart(o)
    if No % oo:
        raise SystemExit('order does not divide group order')
    idx = No // oo
    fo = {q: a for q, a in fa.items() if q != 2}
    EL = VL = 0.0
    hq = 0
    ehq = 0.0
    for q, a in fo.items():
        ge = [q ** -j for j in range(0, a + 1)]   # P(v >= j)
        Ev = sum(ge[1:])
        Ev2 = sum((2 * j - 1) * ge[j] for j in range(1, a + 1))
        EL += Ev * math.log(q)
        VL += (Ev2 - Ev * Ev) * math.log(q) ** 2
        if q <= 50:
            hq += (idx % q == 0)
            ehq += 1.0 / q
    # tail probabilities P(idx >= t): enumerate divisors d of No with d < t
    qs = sorted(fo)

    def pv(q, a, j):
        hi = q ** -j
        lo = q ** -(j + 1) if j < a else 0.0
        return hi - lo
    base = 1.0
    for q in qs:
        base *= pv(q, fo[q], 0)
    pt = []
    tmax = max(TT)
    small = []   # (d, P(idx = d)) for d < tmax

    def rec(i, d, pr):
        if i == len(qs):
            small.append((d, pr))
            return
        q, a = qs[i], fo[qs[i]]
        p0 = pv(q, a, 0)
        dd = d
        for j in range(0, a + 1):
            if dd >= tmax:
                break
            rec(i + 1, dd, pr / p0 * pv(q, a, j) if p0 > 0 else 0.0)
            dd *= q
    rec(0, 1, base)
    for t in TT:
        pt.append(1.0 - sum(pr for d, pr in small if d < t))
    ot = [int(idx >= t) for t in TT]
    U, P = [], []
    for T in TU:
        R = 1
        for q, a in fo.items():
            if q > T:
                R *= q ** a
        U.append(int(all(q <= T for q in fo if oo % q == 0)))
        P.append(1.0 / R)
    return dict(idx=idx, L=math.log(idx), EL=EL, VL=VL, hq=hq, ehq=ehq, ot=ot, pt=pt, U=U, P=P,
                bN=No.bit_length())


def load(files):
    phi, rec, ctl, ends = {}, {}, [], []
    for fn in files:
        for line in open(fn):
            w = line.split()
            if not w:
                continue
            if w[0].startswith('#'):
                ends.append((fn, line.strip()))
                continue
            if w[0] == 'PHI':
                key = (w[1], int(w[2]))
                phi[key] = None if w[3] == 'FAIL' else parse_fac(w[3])
            elif w[0] == 'REC':
                rec[(w[1], int(w[2]), int(w[3]))] = None if w[4] == 'FAILN' else parse_fac(w[4])
            elif w[0] == 'CTL':
                ctl.append((w[1], int(w[2]), int(w[3]), None if w[4] == 'FAILN' else parse_fac(w[4])))
    return phi, rec, ctl, ends


def forced_set(phi, rec, verify=True):
    """re-derive the forced primes from the PHI factorisations; verify every factorisation."""
    out = []
    nfailphi = 0
    for (mode, x), fa in sorted(phi.items()):
        if fa is None:
            nfailphi += 1
            continue
        c = phi2(x) if mode == 'F' else phistar(x)
        if verify:
            check_fac(c, fa, 'Phi %s %d' % (mode, x))
        for p in sorted(fa):
            if p <= 5 or x % p == 0 or jacobi_symbol(5, p) != -1:
                continue
            key = (mode, x, p)
            if key not in rec:
                raise SystemExit('MISSING REC for %s' % (key,))
            out.append(key)
    return out, nfailphi


def gen(argv):
    K, seed, outgp = int(argv[0]), int(argv[1]), argv[2]
    phi, rec, _, _ = load(argv[3:])
    fs, _ = forced_set(phi, rec, verify=False)
    rnd = random.Random(seed)
    items = []
    nskip = 0
    for mode, x, p in fs:
        B = p.bit_length()
        lo, hi = 1 << (B - 1), (1 << B) - 1
        M = 5 * x
        if (hi - lo) // M < 64 * B:
            nskip += 1
            continue
        for _ in range(K):
            r5 = rnd.choice((2, 3))
            rx = 1 % x if mode == 'F' else (-1) % x
            # CRT for residue mod 5x (gcd(x,5) = 1)
            r = next(t for t in range(rx, M, x) if t % 5 == r5)
            kmin = (lo - r + M - 1) // M
            kmax = (hi - r) // M
            for _try in range(100000):
                q = r + M * rnd.randint(kmin, kmax)
                if isprime(q):
                    break
            else:
                raise SystemExit('no control found')
            items.append((1 if mode == 'F' else 2, x, q))
    with open(outgp, 'w') as fh:
        fh.write('CT = [' + ','.join('[%d,%d,%d]' % t for t in items) + '];\n')
    print('forced %d, controls %d, skipped (sparse class) %d' % (len(fs), len(items), nskip))


def summarize(rows):
    n = len(rows)
    if n == 0:
        return None
    s = {'n': n}
    so = sum(r['L'] for r in rows)
    sp = sum(r['EL'] for r in rows)
    sv = sum(r['VL'] for r in rows)
    s['log'] = (so, sp, (so - sp) / math.sqrt(sv) if sv > 0 else float('nan'))
    ho = sum(r['hq'] for r in rows)
    he = sum(r['ehq'] for r in rows)
    s['hq'] = (ho, he, (ho - he) / math.sqrt(he) if he > 0 else float('nan'))
    for k, t in enumerate(TT):
        o = sum(r['ot'][k] for r in rows)
        p = sum(r['pt'][k] for r in rows)
        v = sum(r['pt'][k] * (1 - r['pt'][k]) for r in rows)
        s['t%d' % t] = (o, p, (o - p) / math.sqrt(v) if v > 0 else float('nan'))
    for k, T in enumerate(TU):
        o = sum(r['U'][k] for r in rows)
        p = sum(r['P'][k] for r in rows)
        v = sum(r['P'][k] * (1 - r['P'][k]) for r in rows)
        s['U%d' % T] = (o, p, (o - p) / math.sqrt(v) if v > 0 else float('nan'))
        # non-deterministic part: primes whose group order is NOT T-smooth (P < 1)
        nd = [r for r in rows if r['P'][k] < 1.0]
        s['U%d_nd' % T] = (sum(r['U'][k] for r in nd), sum(r['P'][k] for r in nd), len(nd))
    return s


def stats(argv):
    sets = {}
    files = []
    i = 0
    while i < len(argv):
        if argv[i] == '--set':
            name, fl = argv[i + 1].split('=')
            sets[name] = fl.split(',')
            i += 2
        else:
            files.append(argv[i])
            i += 1
    phi, rec, ctl, ends = load(files)
    for fn, e in ends:
        print('%s: %s' % (fn, e))
    fs, nfailphi = forced_set(phi, rec)
    print('exponents: %d (F %d, L %d), Phi factorisation failures %d'
          % (len(phi), sum(1 for k in phi if k[0] == 'F'), sum(1 for k in phi if k[0] == 'L'), nfailphi))
    rows = []
    nfailN = 0
    usable = []
    for mode, x, p in fs:
        fa = rec[(mode, x, p)]
        if fa is None:
            nfailN += 1
            continue
        if mode == 'F':
            N = p + 1
            check_fac(N, fa, 'p+1 of %d' % p)
            assert pow(2, x, p) == 1 and all(pow(2, x // r, p) != 1 for r in primefactors(x)), 'ord2 != e'
            assert fibmod(N, p) == 0
            o = order_from(N, fa, lambda m: fibmod(m, p) == 0)
            ord2, z = x, o
        else:
            N = p - 1
            check_fac(N, fa, 'p-1 of %d' % p)
            assert fibmod(x, p) == 0 and all(fibmod(x // r, p) != 0 for r in primefactors(x)), 'z != f'
            o = order_from(N, fa, lambda m: pow(2, m, p) == 1)
            ord2, z = o, x
        st = kummer_stats(N, fa, o)
        st.update(tag='FF' if mode == 'F' else 'LF', x=x, p=p, bp=p.bit_length(), ord2=ord2, z=z)
        rows.append(st)
        gx = max(primefactors(oddpart(x)) or [1])
        for k, T in enumerate(TU):
            if st['bp'] > 33 and st['U'][k] and gx <= T:
                usable.append((T, st['tag'], x, p, ord2, z))
    nctl_fail = 0
    for tag, x, p, fa in ctl:
        if fa is None:
            nctl_fail += 1
            continue
        assert isprime(p) and jacobi_symbol(5, p) == -1
        if tag == 'FFc':
            assert p % x == 1 % x
            N = p + 1
            check_fac(N, fa, 'ctl p+1')
            o = order_from(N, fa, lambda m: fibmod(m, p) == 0)
        else:
            assert (p + 1) % x == 0
            N = p - 1
            check_fac(N, fa, 'ctl p-1')
            o = order_from(N, fa, lambda m: pow(2, m, p) == 1)
        st = kummer_stats(N, fa, o)
        st.update(tag=tag, x=x, p=p, bp=p.bit_length())
        rows.append(st)
    print('forced primes: FF %d, LF %d; other-order factoring failures %d; controls %d (fail %d)'
          % (sum(r['tag'] == 'FF' for r in rows), sum(r['tag'] == 'LF' for r in rows), nfailN,
             sum(r['tag'] in ('FFc', 'LFc') for r in rows), nctl_fail))
    for tg in ('FF', 'LF'):
        b = sorted(r['bp'] for r in rows if r['tag'] == tg)
        if b:
            print('  %s bits: min %d median %d max %d; n(p>2^33) %d, n(p>2^64) %d'
                  % (tg, b[0], b[len(b) // 2], b[-1], sum(x > 33 for x in b), sum(x > 64 for x in b)))
    strata = [('all', lambda r: True), ('p<2^33', lambda r: r['bp'] <= 33), ('p>2^33', lambda r: r['bp'] > 33),
              ('p>2^64', lambda r: r['bp'] > 64)]
    keys = ['log', 'hq'] + ['t%d' % t for t in TT] + ['U%d' % T for T in TU]
    for side, cside in (('FF', 'FFc'), ('LF', 'LFc')):
        for sname, sel in strata:
            a = summarize([r for r in rows if r['tag'] == side and sel(r)])
            c = summarize([r for r in rows if r['tag'] == cside and sel(r)])
            if a is None:
                continue
            print('%s %-7s n=%d  ctl n=%d' % (side, sname, a['n'], c['n'] if c else 0))
            for k in keys:
                o, p, z = a[k]
                line = '   %-6s obs %8.2f pred %8.2f z %+6.2f ratio %6.3f' % (k, o, p, z, o / p if p else float('nan'))
                if c:
                    o2, p2, z2 = c[k]
                    r2 = o2 / p2 if p2 else float('nan')
                    line += ' | ctl obs %8.2f pred %8.2f z %+6.2f ratio %6.3f' % (o2, p2, z2, r2)
                    if r2 and r2 == r2 and p:
                        line += ' kappa %6.3f' % ((o / p) / r2)
                if k.startswith('U'):
                    o3, p3, n3 = a[k + '_nd']
                    line += ' | nondet: n %d obs %d pred %.3g' % (n3, o3, p3)
                print(line)
    print('usable forced primes with p > 2^33 (both odd orders T-smooth):')
    for t in usable:
        T, tg, x, p, o2, z = t
        print('  T=%d %s x=%d p=%d (%d bits, p mod 8 = %d) ord2=%d z=%d class (v2 ord2, v2 z)=(%d,%d) '
              'gpf odd ord2 %d, gpf odd z %d' % (T, tg, x, p, p.bit_length(), p % 8, o2, z, v2(o2), v2(z),
                                                 max(primefactors(oddpart(o2)) or [1]),
                                                 max(primefactors(oddpart(z)) or [1])))
    # dump per-prime forced records for later cross-checks
    with open('adv_forced_records.txt', 'w') as fh:
        for r in rows:
            if r['tag'] in ('FF', 'LF'):
                fh.write('%s %d %d %d %d %d %d %.6g %.6g %.6g %d %.6g %s %s %s %s\n'
                         % (r['tag'], r['x'], r['p'], r['bp'], r['ord2'], r['z'], r['idx'], r['L'], r['EL'], r['VL'],
                            r['hq'], r['ehq'], ','.join(map(str, r['ot'])), ','.join('%.6g' % v for v in r['pt']),
                            ','.join(map(str, r['U'])), ','.join('%.6g' % v for v in r['P'])))


if __name__ == '__main__':
    if sys.argv[1] == 'gen':
        gen(sys.argv[2:])
    else:
        stats(sys.argv[2:])
