#!/usr/bin/env python3
"""xn line: exact (Fraction) check at k = 9 (random rational members, with/without 2t+1) of
 (2) V''(-1) = (k-1)(k-2) 2^(1-k) + 2^(-k) e2(beta),   (2r) same for the reflected family with gamma,
 (7) even-subset identity V(t) = 2^(1-k) sum_j t^(k-1-2j) (t+1)^(2j) e_2j(h_i(t)/g_i(t)) at rational t,
members sampled inside the real domain (so V >= 2^(1-k) t^(k-1) must hold), where V := (F-1)/(t P) as a rational function (Taylor data from xn_core.G_tau: (1-tau)^(k-1) G = V exactly)."""
import random
from fractions import Fraction as Fr
from itertools import combinations
from xn_core import G_tau

random.seed(9)


def e(vals, j):
    s = Fr(0)
    for c in combinations(vals, j):
        p = Fr(1)
        for x in c:
            p *= x
        s += p
    return s


for lin in (0, 1):
    k, m = 9, 9 - lin
    for trial in range(2):
        A, C = [], []
        while len(A) < m:      # sample INSIDE the real domain c^2 < 8a, (a - c)^2 < 8a (g, h without real roots)
            a = Fr(random.randint(1, 3100), 100)
            c = Fr(random.randint(-1600, 1600), 100)
            if c * c < 8 * a and (a - c) ** 2 < 8 * a:
                A.append(a); C.append(c)
        for refl in (0, 1):
            CC = [a - c for a, c in zip(A, C)] if refl else C
            V2 = G_tau(A, CC, lin, k, 3)[2]                  # [tau^2] V = V''(-1)/2
            beta = [2 + a - c for a, c in zip(A, CC)] + [Fr(2)] * lin
            ok2 = 2 * V2 == Fr((k - 1) * (k - 2), 2 ** (k - 1)) + e(beta, 2) / 2 ** k
            print(f"lin={lin} trial {trial} refl={refl}: (2) exact: {ok2}")
        # (7) at a few rational t
        for t in (Fr(-7, 2), Fr(-3, 2), Fr(-3, 4), Fr(-1, 3), Fr(1, 3), Fr(5, 2)):
            f = [2 * t + 1 + t * (t + 1) * (a * t + c) for a, c in zip(A, C)] + [2 * t + 1] * lin
            F = Fr(1)
            for x in f:
                F *= x
            g = [(x - 1) / t for x in f]
            h = [(x + 1) / (t + 1) for x in f]
            P = Fr(1)
            for x in g:
                P *= x
            V = (F - 1) / (t * P)
            r = [hh / gg for hh, gg in zip(h, g)]
            rhs = sum(t ** (k - 1 - 2 * j) * (t + 1) ** (2 * j) * e(r, 2 * j) for j in range((k - 1) // 2 + 1)) / 2 ** (k - 1)
            print(f"   (7) at t = {t}: exact {V == rhs};  V >= 2^-8 t^8: {V >= t ** 8 / 256}")
