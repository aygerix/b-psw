#!/usr/bin/env python3
"""xn line: EXACT cross-check (Fractions) of the k = 9 moment forms used by xn_lp9.py / xn_farkas9.py against the
direct series computation of xn_core (Vpoly from the Laurent expansion).  Random rational configurations."""
from fractions import Fraction as Fr
import random
from xn_core import vpoly_s
import xn_farkas9 as X

random.seed(5)


def moments(U, Vv):
    w = [u - Fr(1, 2) for u in U]
    M = lambda p, j: sum(v ** p * x ** j for v, x in zip(Vv, w))
    return dict(S=M(1, 0), A1=M(1, 1), A2=M(1, 2), A4=M(1, 4), B0=M(2, 0), B2=M(2, 2), C0=M(3, 0))


def kval(K, mo):
    S = mo["S"]
    val = K.get("SS", 0) * S * S + K.get("S", 0) * S + K.get("c", 0) + K.get("SA2", 0) * S * mo["A2"] + K.get("SB0", 0) * S * mo["B0"]
    val += K.get("P", 0) * mo["A1"] ** 2 / 2
    for key in ("A2", "B0", "A4", "B2", "C0"):
        val += K.get(key, 0) * mo[key]
    return val


for lin in (0, 1):
    m, K0, K1, KP, TE, T12 = X.model(lin)
    k = 9
    for trial in range(3):
        U = [Fr(random.randint(-30, 30), random.randint(1, 9)) for _ in range(m)]
        Vv = [Fr(random.randint(1, 40), random.randint(1, 9)) for _ in range(m)]
        T = vpoly_s(U, Vv, lin, k)
        Tr = vpoly_s([1 - u for u in U], Vv, lin, k)
        # V(t) = sum_n T_n t^(8-n)
        V = lambda TT, t: sum(TT[n] * t ** (8 - n) for n in range(9))
        dV = lambda TT, t: sum(TT[n] * (8 - n) * t ** (7 - n) for n in range(8))
        d2V = lambda TT, t: sum(TT[n] * (8 - n) * (7 - n) * t ** (6 - n) for n in range(7))
        E0 = (V(T, -1) - Fr(1, 256) + V(Tr, -1) - Fr(1, 256)) / 2
        E1 = (dV(T, -1) + Fr(8, 256) + dV(Tr, -1) + Fr(8, 256)) / 2
        PH = (d2V(T, -1) - Fr(56, 256) + d2V(Tr, -1) - Fr(56, 256)) / 2
        mo = moments(U, Vv)
        ok = [E0 == kval(K0, mo), E1 == kval(K1, mo), PH == kval(KP, mo)]
        okT = [(T[n] + Tr[n]) / 2 == kval(TE[n - 3], mo) for n in range(3, 9)]
        okT12 = [(T[1] + Tr[1]) / 2 == T12[0], (T[2] + Tr[2]) / 2 == T12[1]]
        print(f"lin={lin} trial {trial}: Es0, Es1, Phi exact match: {ok}; T3e..T8e: {okT}; T1e,T2e: {okT12}")
