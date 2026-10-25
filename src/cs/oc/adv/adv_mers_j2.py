#!/usr/bin/env python3
"""Adversarial check (line OC): Mersenne Lucas classes, Theorem M, class collisions.
Independent of oc_j2table.gp: arithmetic in F_p[r]/(r^2 - D) (r = sqrt D), gamma = alpha/beta = alpha^2/Q with
alpha = (1 + r)/2; j2 = v2(ord gamma) found by repeated squaring of gamma^(odd part of group order).
Usage: adv_mers_j2.py QMAX [table_to_compare]"""
import sys, time
from sympy import jacobi_symbol

MEXP = [3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423]

def mul(a, b, D, p):
    # (a0 + a1 r)(b0 + b1 r), r^2 = D
    return ((a[0]*b[0] + D*a[1]*b[1]) % p, (a[0]*b[1] + a[1]*b[0]) % p)

def pw(a, e, D, p):
    r = (1, 0)
    while e:
        if e & 1: r = mul(r, a, D, p)
        a = mul(a, a, D, p); e >>= 1
    return r

def v2(x): return (x & -x).bit_length() - 1

def j2_of(p, D):
    Q = (1 - D) // 4
    eps = int(jacobi_symbol(D % p, p))
    if eps == 0: return 0, None
    if Q % p == 0: return eps, -1          # degenerate: p | Q (U_k = 1 for all k >= 1), never a BPSW factor
    i2 = pow(2, -1, p)
    alpha = (i2, i2)
    gam = mul(alpha, alpha, D, p)
    iq = pow(Q % p, -1, p)
    gam = (gam[0] * iq % p, gam[1] * iq % p)
    N = p - eps
    m = N >> v2(N)
    g = pw(gam, m, D, p)
    j = 0
    while g != (1, 0):
        g = mul(g, g, D, p); j += 1
        if j > v2(N): raise RuntimeError('order does not divide group order')
    return eps, j

def selfridge_nonsquare(limit):
    out, D = [], 5
    while abs(D) <= limit:
        if not (D > 0 and int(D**0.5 + 0.5)**2 == D): out.append(D)
        D = -D - 2 if D > 0 else -D + 2
    return out

if __name__ == '__main__':
    QMAX = int(sys.argv[1]); t0 = time.time()
    Ds = selfridge_nonsquare(103)
    res = {}
    for q in MEXP:
        if q > QMAX: break
        p = (1 << q) - 1
        for D in Ds:
            e, j = j2_of(p, D)
            res[(D, q)] = (e, j)
        print(f'q={q} done {time.time()-t0:.1f}s', flush=True)
    # compare with table
    if len(sys.argv) > 2:
        nm = nc = 0
        for line in open(sys.argv[2]):
            D, q, e, j = map(int, line.split())
            if (D, q) not in res: continue
            nc += 1
            e2, j2 = res[(D, q)]
            if e2 == 0: e2, j2 = 0, -9
            if (e, j) != (e2, j2): nm += 1; print('TABLE MISMATCH', D, q, (e, j), (e2, j2))
        print(f'compared {nc} table rows, mismatches {nm}')
    # Theorem M
    bad = 0; cnt = {}
    for (D, q), (e, j) in res.items():
        if e != -1 or j == -1: continue
        if D in (5, -11, 17):
            cnt[D] = cnt.get(D, 0) + 1
            if j != q: bad += 1; print('THM M FAIL', D, q, j)
        if D in (-7, -15):
            cnt[D] = cnt.get(D, 0) + 1
            if j != q - 1: bad += 1; print('THM M FAIL', D, q, j)
        if D == -7 and q % 6 != 5: bad += 1; print('D=-7 inert with q !=5 mod 6', q)
        if D == -15 and q % 4 != 3: bad += 1; print('D=-15 inert with q !=3 mod 4', q)
    print('Theorem M inert cases checked per D:', cnt, 'failures', bad)
    # D = -31 closed form, OC4 law on every inert (Q/p)=1 case
    noc4 = bad4 = n31 = bad31 = 0
    for (D, q), (e, j) in res.items():
        if e != -1 or j == -1: continue
        p = (1 << q) - 1; Q = (1 - D) // 4
        if jacobi_symbol(Q % p, p) != 1:
            if j != q: print('Lemma 18 FAIL', D, q, j)
            continue
        s = pow(Q % p, (p + 1) // 4, p); assert s * s % p == Q % p
        sym = jacobi_symbol((2*Q + s) % p, p); sym2 = jacobi_symbol((2*Q - s) % p, p)
        noc4 += 1
        if ((j == q - 1) != (sym == -1)) or sym != sym2: bad4 += 1; print('OC4 FAIL', D, q, j, sym, sym2)
        if D == -31 and q >= 9:
            n31 += 1; m = (q - 5) // 2
            pred = jacobi_symbol(pow(2, m, 31) + 1, 31)
            if (pred == -1) != (j == q - 1): bad31 += 1; print('D=-31 closed form FAIL', q, j, pred)
    print(f'OC4 on Mersenne: {noc4} inert (Q/p)=1 cases, failures {bad4}; D=-31 closed form {n31} cases, failures {bad31}')
    # collisions among inert members, per D
    for D in Ds:
        js = {}
        for q in MEXP:
            if (D, q) in res and res[(D, q)][0] == -1 and res[(D, q)][1] != -1: js.setdefault(res[(D, q)][1], []).append(q)
        for j, qs in js.items():
            if len(qs) > 1: print('INERT CLASS COLLISION D', D, 'j2', j, 'q', qs)
    # split members: j2 <= 1 always
    print('split j2 max:', max(j for (D, q), (e, j) in res.items() if e == 1))
    print(f'total {time.time()-t0:.1f}s')
