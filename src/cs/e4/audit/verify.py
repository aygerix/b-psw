#!/usr/bin/env python3
"""Independent audit of the E4 (FIC) records.

Reads the instrumented rerun outputs <run>.out / <run>.rec / <run>.phi produced by fic_aud.gp, and
  (1) checks the .out records are byte-identical to the original out/ files (done outside, by cmp);
  (2) re-derives every forced prime from an independently computed Phi_e(2) / Phi*_f and the logged
      factorisation (product check + sympy.isprime on every factor), and checks completeness of the
      exponent loop;
  (3) re-verifies forcing (ord_p 2 = e, z(p) = f exactly), every control's class / size / inertness,
      every logged factorisation of the other group order, and recomputes the other order with code
      written here (own Fibonacci doubling, own order routine);
  (4) recomputes every per-prime statistic with independently coded Kummer formulas and compares with
      the record fields;
  (5) recomputes the summary statistics (obs, pred, z) per tag and stratum from the recomputed values.
Usage: verify.py run1 [run2 ...]   (run = basename, e.g. aF)
"""
import math
import sys
from sympy import isprime, jacobi_symbol, factorint, divisors

TT = [3, 10, 100, 1000]
TU = [200, 1500, 10000]


def parse_fac(s):
    if s == '1':
        return {}
    out = {}
    for t in s.split('*'):
        q, a = t.split('^')
        out[int(q)] = out.get(int(q), 0) + int(a)
    return out


def prod_fac(fa):
    r = 1
    for q, a in fa.items():
        r *= q ** a
    return r


def mobius(n):
    f = factorint(n)
    if any(a > 1 for a in f.values()):
        return 0
    return -1 if len(f) % 2 else 1


def eulerphi(n):
    r = n
    for q in factorint(n):
        r = r // q * (q - 1)
    return r


def fib_pair(n, m):
    """(F_n, F_{n+1}) mod m by fast doubling."""
    if n == 0:
        return (0, 1 % m)
    a, b = fib_pair(n >> 1, m)
    c = (a * ((2 * b - a) % m)) % m
    d = (a * a + b * b) % m
    if n & 1:
        return (d, (c + d) % m)
    return (c, d)


def fibmod(n, m):
    return fib_pair(n, m)[0]


def fib(n):
    return fib_pair(n, 1 << (n + 2))[0] if n > 0 else 0


def phi_e_2(e):
    num, den = 1, 1
    for d in divisors(e):
        mu = mobius(e // d)
        if mu == 1:
            num *= (1 << d) - 1
        elif mu == -1:
            den *= (1 << d) - 1
    assert num % den == 0
    return num // den


def phistar_f(f):
    num, den = 1, 1
    for d in divisors(f):
        mu = mobius(f // d)
        if mu == 1:
            num *= fib(d)
        elif mu == -1:
            den *= fib(d)
    assert num % den == 0
    return num // den


def order_generic(N, fa, is_one):
    """smallest n | N with is_one(n), given factorisation fa of N (is_one is closed under multiples)."""
    n = N
    for q, a in fa.items():
        for _ in range(a):
            if is_one(n // q):
                n //= q
            else:
                break
    return n


def ord2(p, fa):
    return order_generic(p - 1, fa, lambda n: pow(2, n, p) == 1)


def zent(p, fa):
    return order_generic(p + 1, fa, lambda n: fibmod(n, p) == 0)


def oddpart(n):
    while n % 2 == 0:
        n //= 2
    return n


def stats(N, fa, o):
    """independent recomputation of the record statistics"""
    No, oo = oddpart(N), oddpart(o)
    assert No % oo == 0
    idx = No // oo
    oq = {q: a for q, a in fa.items() if q != 2}
    Elog = Vlog = 0.0
    for q, a in oq.items():
        Ev = sum(q ** -j for j in range(1, a + 1))           # sum_j P(v >= j)
        Ev2 = sum((2 * j - 1) * q ** -j for j in range(1, a + 1))
        Elog += Ev * math.log(q)
        Vlog += (Ev2 - Ev * Ev) * math.log(q) ** 2
    hq = sum(1 for q in oq if q <= 50 and idx % q == 0)
    ehq = sum(1.0 / q for q in oq if q <= 50)

    def pv(q, a, j):  # P(v_q(index) = j)
        return q ** -j * (1 - 1.0 / q) if j < a else float(q) ** -a

    ot, pt = [], []
    for t in TT:
        # enumerate odd divisors d < t of No and sum P(index = d) = prod_q P(v_q = v_q(d))
        s = 0.0
        stack = [(1, 0)]
        qs = sorted(oq)
        # recursive enumeration of divisors below t
        def rec(i, d, pr):
            nonlocal s
            if i == len(qs):
                s += pr
                return
            q, a = qs[i], oq[qs[i]]
            dd = d
            for j in range(0, a + 1):
                if dd >= t:
                    break
                rec(i + 1, dd, pr * pv(q, a, j))
                dd *= q
        rec(0, 1, 1.0)
        pt.append(1 - s)
        ot.append(1 if idx >= t else 0)
    U, P = [], []
    for T in TU:
        sT = 1
        for q, a in oq.items():
            if q <= T:
                sT *= q ** a
        U.append(1 if all(q <= T for q in factorint(oo)) else 0)
        P.append(sT / No)
    return dict(bN=No.bit_length(), L=math.log(idx), EL=Elog, VL=Vlog, hq=hq, ehq=ehq, o=ot, p=pt, U=U, P=P, idx=idx)


def close(a, b, rel=2e-5, absol=1e-5):
    return abs(a - b) <= max(absol, rel * max(abs(a), abs(b)))


def load_out(path):
    rows = []
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        f = line.split()
        rows.append(f)
    return rows


def main():
    runs = sys.argv[1:]
    allrec = []
    nerr = 0
    for run in runs:
        outrows = load_out(run + '.out')
        recs = [l.split() for l in open(run + '.rec') if l.strip()]
        phis = [l.split() for l in open(run + '.phi') if l.strip()]
        assert len(outrows) == len(recs), (run, len(outrows), len(recs))
        # ---- (2) forced sets from the logged full factorisation of Phi
        forced_from_phi = {}
        exps_done = {'F': set(), 'L': set()}
        fails = {'F': set(), 'L': set()}
        for ph in phis:
            if ph[0] in ('PF', 'PL'):
                side = ph[0][1]
                ef = int(ph[1])
                fa = parse_fac(ph[2])
                C = phi_e_2(ef) if side == 'F' else phistar_f(ef)
                if prod_fac(fa) != C:
                    print('ERR product of Phi factorisation', run, ph[:2]); nerr += 1
                for q in fa:
                    if not isprime(q):
                        print('ERR composite factor in Phi', run, ph[:2], q); nerr += 1
                sel = sorted(q for q in fa if q > 5 and ef % q != 0 and jacobi_symbol(5, q) == -1)
                forced_from_phi[(side, ef)] = sel
                exps_done[side].add(ef)
            elif ph[0] in ('FAILF', 'FAILL'):
                fails[ph[0][-1]].add(int(ph[1]))
            elif ph[0] in ('FAILNF', 'FAILNL'):
                print('NOTE FAILN', run, ph)
        forced_from_rec = {}
        last_forced = None
        for k, (r, o) in enumerate(zip(recs, outrows)):
            tag, ef, p, ordr = r[0], int(r[1]), int(r[2]), int(r[3])
            fa = parse_fac(r[4])
            if o[0] != tag or int(o[1]) != ef or int(o[2]) != p.bit_length():
                print('ERR rec/out mismatch', run, k, r[:3], o[:3]); nerr += 1
            if not isprime(p):
                print('ERR p not prime', run, r[:3]); nerr += 1
            if jacobi_symbol(5, p) != -1:
                print('ERR p not inert', run, r[:3]); nerr += 1
            N = p + 1 if tag.startswith('FF') else p - 1
            if prod_fac(fa) != N or not all(isprime(q) for q in fa):
                print('ERR bad factorisation of other group order', run, r[:3]); nerr += 1
            if tag == 'FF':
                # forcing: ord_p 2 = e exactly
                ok = pow(2, ef, p) == 1 and all(pow(2, ef // q, p) != 1 for q in factorint(ef))
                if not ok:
                    print('ERR ord_p 2 != e', run, r[:3]); nerr += 1
                forced_from_rec.setdefault(('F', ef), []).append(p)
                last_forced = (tag, ef, p)
                oth = zent(p, fa)
                forcedside = ef
            elif tag == 'LF':
                ok = fibmod(ef, p) == 0 and all(fibmod(ef // q, p) != 0 for q in factorint(ef))
                if not ok:
                    print('ERR z(p) != f', run, r[:3]); nerr += 1
                forced_from_rec.setdefault(('L', ef), []).append(p)
                last_forced = (tag, ef, p)
                oth = ord2(p, fa)
                forcedside = ef
            elif tag == 'FFc':
                if not (last_forced and last_forced[0] == 'FF' and last_forced[1] == ef
                        and last_forced[2].bit_length() == p.bit_length() and p % ef == 1):
                    print('ERR FF control class/size', run, r[:3], last_forced); nerr += 1
                oth = zent(p, fa)
                forcedside = None
            elif tag == 'LFc':
                if not (last_forced and last_forced[0] == 'LF' and last_forced[1] == ef
                        and last_forced[2].bit_length() == p.bit_length() and (p + 1) % ef == 0):
                    print('ERR LF control class/size', run, r[:3], last_forced); nerr += 1
                oth = ord2(p, fa)
                forcedside = None
            if oth != ordr:
                print('ERR other-side order differs', run, r[:3], oth, ordr); nerr += 1
            st = stats(N, fa, oth)
            # compare with record fields
            vals = [float(x) for x in o[3:]]
            mine = [st['bN'], st['L'], st['EL'], st['VL'], st['hq'], st['ehq']] + st['o'] + st['p'] + \
                   [st['U'][0], st['P'][0], st['U'][1], st['P'][1], st['U'][2], st['P'][2]]
            for i, (a, b) in enumerate(zip(vals, mine)):
                if not close(a, b):
                    print('ERR stat field %d' % i, run, r[:3], a, b); nerr += 1
                    break
            st.update(tag=tag, ef=ef, p=p, bp=p.bit_length(), ordr=oth, run=run, N=N, fa=fa)
            allrec.append(st)
        # completeness: forced primes from Phi factorisations == forced primes in records
        for key, sel in forced_from_phi.items():
            got = sorted(forced_from_rec.get(key, []))
            if got != sel:
                print('ERR forced set mismatch', run, key, sel, got); nerr += 1
        for key in forced_from_rec:
            if key not in forced_from_phi:
                print('ERR forced prime with no Phi line', run, key); nerr += 1
        # completeness of the exponent loop, from the driver parameters
        import re
        drv = open('drv_' + run + '.gp').read()
        par = dict(EMIN=2, FMIN=3, PHIBITS=10 ** 9, PHILO=0)
        for k, v in re.findall(r'(\w+)=(\d+);', drv):
            par[k] = int(v)
        for e in range(par['EMIN'], par['EMAX'] + 1):
            if e % 5 == 0 or not (par['PHILO'] < eulerphi(e) <= par['PHIBITS']):
                continue
            if e not in exps_done['F'] and e not in fails['F']:
                print('ERR exponent e not processed', run, e); nerr += 1
        for f in range(par['FMIN'], par['FMAX'] + 1):
            if f % 5 == 0 or not (par['PHILO'] < 0.6942 * eulerphi(f) <= par['PHIBITS']):
                continue
            if f not in exps_done['L'] and f not in fails['L']:
                print('ERR exponent f not processed', run, f); nerr += 1
        n_extra = len([e for e in exps_done['F'] if not (par['PHILO'] < eulerphi(e) <= par['PHIBITS'])]) + \
                  len([f for f in exps_done['L'] if not (par['PHILO'] < 0.6942 * eulerphi(f) <= par['PHIBITS'])])
        if n_extra:
            print('ERR exponents outside band processed', run, n_extra); nerr += 1
        print('   params', par)
        print('run %s: %d records, %d Phi lines (F %d, L %d), Phi failures F %d L %d'
              % (run, len(recs), len(forced_from_phi), len(exps_done['F']), len(exps_done['L']),
                 len(fails['F']), len(fails['L'])))
        print('   failed exponents:', sorted(fails['F']), sorted(fails['L']))
    print('TOTAL ERRORS:', nerr)
    return allrec


if __name__ == '__main__':
    import pickle
    recs = main()
    pickle.dump(recs, open('verified_' + '_'.join(sys.argv[1:]) + '.pkl', 'wb'))
