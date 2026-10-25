#!/usr/bin/env python3
"""aud13: exact checks of the two identities the cuts rest on (k = 9, random rational members, both cases):
 (P) H(t) := (F - 1)/(t prod g) = 2^(1-k) sum_j t^(k-1-2j) (t+1)^(2j) e_(2j)(h_i/g_i)   at random rational t
     (gives POS: H(t) >= 2^-8 t^8 on R when all g_i, h_i > 0);
 (D) H(-1) = 2^-8, H'(-1) = -8 2^-8, H''(-1) = 2^-8 [56 + e_2(beta)/2], beta_i = h_i(-1)   (gives PHI and the N3 relation).
Exact rational arithmetic; no moment forms.
"""
import random
from fractions import Fraction as Fr


def pmul(a, b):
    o = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            o[i + j] += x * y
    return o


def pev(p, x):
    r = Fr(0)
    for c in reversed(p):
        r = r * x + c
    return r


def pder(p):
    return [i * p[i] for i in range(1, len(p))] or [Fr(0)]


def esym(xs, n):
    e = [Fr(1)] + [Fr(0)] * n
    for x in xs:
        for j in range(n, 0, -1):
            e[j] += e[j - 1] * x
    return e


rnd = random.Random(5)
K = 9
bad = 0
for trial in range(40):
    lin = trial % 2
    m = K - lin
    mem = [(Fr(rnd.randint(1, 300), rnd.randint(1, 60)), Fr(rnd.randint(-200, 200), rnd.randint(1, 50))) for _ in range(m)]
    F, D = [Fr(1)], [Fr(0), Fr(1)]
    gs, hs = [], []
    for a, c in mem:
        F = pmul(F, [Fr(1), c + 2, a + c, a])
        g = [c + 2, a + c, a]
        h = [Fr(2), c, a]
        D = pmul(D, g)
        gs.append(g)
        hs.append(h)
    if lin:
        F = pmul(F, [Fr(1), Fr(2)])
        D = pmul(D, [Fr(2)])
        gs.append([Fr(2)])
        hs.append([Fr(2)])
    N = list(F)
    N[0] -= 1
    # (P)
    for _ in range(3):
        t = Fr(rnd.randint(-50, 50), rnd.randint(1, 13))
        if t == 0 or t == -1 or pev(D, t) == 0:
            continue
        H = pev(N, t) / pev(D, t)
        r = [pev(h, t) / pev(g, t) for g, h in zip(gs, hs)]
        e = esym(r, K - 1)
        rhs = Fr(2, 2 ** K) * sum(t ** (K - 1 - 2 * j) * (t + 1) ** (2 * j) * e[2 * j] for j in range((K - 1) // 2 + 1))
        if H != rhs:
            bad += 1
            print("identity (P) FAILS", trial, t)
    # (D): derivatives of H = N/D at -1
    x = Fr(-1)
    N0, N1, N2 = pev(N, x), pev(pder(N), x), pev(pder(pder(N)), x)
    D0, D1, D2 = pev(D, x), pev(pder(D), x), pev(pder(pder(D)), x)
    H0 = N0 / D0
    H1 = (N1 * D0 - N0 * D1) / D0 ** 2
    H2 = (N2 * D0 ** 2 - 2 * N1 * D0 * D1 - N0 * D0 * D2 + 2 * N0 * D1 ** 2) / D0 ** 3
    beta = [pev(h, x) for h in hs]
    e2 = esym(beta, 2)[2]
    ok = (H0 == Fr(1, 256), H1 == Fr(-8, 256), H2 == Fr(1, 256) * (56 + e2 / 2))
    gm1 = [pev(g, x) for g in gs]
    if not all(ok) or any(v != 2 for v in gm1):
        bad += 1
        print("identity (D) FAILS", trial, ok, gm1)
print(f"aud13_hid: 40 random member sets (20 per case): {'ALL EXACT IDENTITIES HOLD' if bad == 0 else str(bad) + ' FAILURES'}")
