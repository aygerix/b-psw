#!/usr/bin/env python3
"""aud13_A_cuts.py -- audit of candidate Theorem 13, items (ii) and (iv). New code (system python3, Fractions only).

(1) Coverage/containment bookkeeping of the certificate files xn_cert9_l{0,1}.gp (read as data):
    contiguity of the R intervals from R = 0, the S range S(R) = (R + 1/2)^2/8 + (m-1)/32, the row kinds and
    parameters used, z >= 0, and that the S-bound (aud13_A_mom.out) is below the cover end.
(2) Direction of every inequality row, tested on the H side (what a family satisfies, since V = H for a family):
    at random exact-rational member sets in the open domain D (u^2 < 8v, (1-u)^2 < 8v), including stress sets
    near the corner v = 1/32 (a -> 32, b -> 18) and with large v:
      CS   (sum v w)^2 <= (sum v)(sum v w^2)
      PHI  Phi_H := 2^-10 (e2(beta) + e2(gamma)) >= 0 and Phi_H == He''(-1) - 56 2^-8 (exact Taylor, finite differences
           are NOT used: Phi_H is compared with the second-order formula only through the identity file; here we check
           the cut algebra)
      N3L  2^10 Phi_H >= 2B sum b - B^2 - sum b^2 - sum om^2 for every B used by a certificate and random B
      N3U  2^10 Phi_H <= Bbar sum b - sum b^2 + (Mtot - 1) sum om^2, Bbar = 18 m + 2 lin
      POS  He(t0) - 2^-8 t0^8 >= 0 at every t0 used by a certificate, He = (H + H~)/2, H(t0) = (F(t0)-1)/(t0 prod g(t0))
    plus a negative control (a deliberately wrong N3U with Bbar = 6 m) that must fail somewhere.
"""
import random
import re
import sys
from fractions import Fraction as F
from math import comb

random.seed(20261025)


def load_certs(fn):
    out = []
    for line in open(fn):
        line = line.strip()
        if not line:
            continue
        txt = re.sub(r"(-?\d+)/(\d+)", r"F(\1,\2)", line)
        out.append(eval(txt, {"__builtins__": {}}, {"F": F}))
    return out


def SofR(R, m):
    return (R + F(1, 2)) ** 2 / 8 + F(m - 1, 32)


def coverage(fn, lin, Sbound):
    m = 9 - lin
    C = load_certs(fn)
    kinds, ok = {}, True
    C2 = sorted(C, key=lambda c: c[1])
    reach = F(0)
    for c in C2:
        assert c[0] == lin, "lin mismatch"
        if c[1] != reach:
            print(f"  GAP/OVERLAP at R = {reach} vs next Ra = {c[1]}")
            ok = False
        if c[2] <= c[1]:
            print("  empty interval", c[1], c[2]); ok = False
        reach = max(reach, c[2])
        if len(c[3]) != 4:
            print("  y length", len(c[3])); ok = False
        for (kd, prm, z) in c[4]:
            if z < 0:
                print("  NEGATIVE z", kd, prm, z); ok = False
            kinds.setdefault(kd, set()).add(prm)
    Send = SofR(reach, m)
    print(f"{fn}: {len(C)} certificates, contiguous from R = 0: {ok and C2[0][1] == 0}; R reach = {reach}; "
          f"S covered = [{SofR(F(0), m)}, {Send}] = [{float(SofR(F(0), m)):.6f}, {float(Send):.6f}]")
    print(f"  S lower end equals m/32 = {F(m, 32)}: {SofR(F(0), m) == F(m, 32)}; S-bound {Sbound:.7f} < cover end: {Sbound < float(Send)}")
    print("  row kinds used (count of distinct parameters):", {k: len(v) for k, v in kinds.items()})
    return kinds


def e2(L):
    s1 = sum(L)
    return (s1 * s1 - sum(x * x for x in L)) / 2


def Hval(mem, t0):
    """H(t0) = (F(t0) - 1)/(t0 prod g(t0)), mem = list of (a, c); lin member is (None, None)"""
    Fv, G = F(1), F(1)
    for (a, c) in mem:
        if a is None:
            f, g = 2 * t0 + 1, F(2)
        else:
            f = 2 * t0 + 1 + t0 * (t0 + 1) * (a * t0 + c)
            g = a * t0 * t0 + (a + c) * t0 + c + 2
        Fv *= f
        G *= g
    if t0 == 0:   # limit: sum_i prod_{j != i} 1/g_j(0)
        gs = [F(2) if a is None else c + 2 for (a, c) in mem]
        return sum(F(1) / prod_excl(gs, i) for i in range(len(gs)))
    return (Fv - 1) / (t0 * G)


def prod_excl(gs, i):
    p = F(1)
    for j, g in enumerate(gs):
        if j != i:
            p *= g
    return p


def rand_member(mode):
    """exact rational (u, v) in the open domain D"""
    while True:
        if mode == "corner":
            v = F(1, 32) + F(random.randint(1, 1000), 10 ** 6)
        elif mode == "big":
            v = F(random.randint(1, 4000), 100)
        else:
            v = F(1, 32) + F(random.randint(1, 3000), 1000)
        u = F(random.randint(-3000, 4000), 1000)
        if u * u < 8 * v and (1 - u) ** 2 < 8 * v:
            return u, v


def cut_tests(lin, kinds, ntr=300):
    m = 9 - lin
    Mtot = m + lin
    Bbar = 18 * m + 2 * lin
    Bs = sorted(kinds.get("N3L", set()))
    t0s = sorted(kinds.get("POS", set()))
    viol = {"CS": 0, "PHI": 0, "N3L": 0, "N3U": 0, "POS": 0, "NEGCTL": 0}
    worst = {"N3U": None, "POS": None}
    for tr in range(ntr):
        mode = ["corner", "big", "mid"][tr % 3]
        us_vs = [rand_member(mode if random.random() < 0.7 else "mid") for _ in range(m)]
        mem = [(1 / v, u / v) for (u, v) in us_vs]
        memr = [(a, a - c) for (a, c) in mem]
        if lin:
            mem = mem + [(None, None)]
            memr = memr + [(None, None)]
        beta = [F(2) if a is None else 2 + a - c for (a, c) in mem]
        gamma = [F(2) if a is None else 2 + c for (a, c) in mem]
        assert all(x > 0 for x in beta + gamma) and all(x < 18 for x in beta + gamma)
        PhiH = (e2(beta) + e2(gamma)) / 1024
        b = [F(2) if a is None else 2 + a / 2 for (a, c) in mem]
        om = [F(0) if a is None else c - a / 2 for (a, c) in mem]
        sb, sb2, so2 = sum(b), sum(x * x for x in b), sum(x * x for x in om)
        assert 1024 * PhiH == sb * sb + sum(om) ** 2 - sb2 - so2
        S = sum(v for (u, v) in us_vs)
        M11 = sum(v * (u - F(1, 2)) for (u, v) in us_vs)
        A2 = sum(v * (u - F(1, 2)) ** 2 for (u, v) in us_vs)
        if not M11 * M11 <= S * A2:
            viol["CS"] += 1
        if not PhiH >= 0:
            viol["PHI"] += 1
        for B in Bs + [F(random.randint(-500, 500), 7) for _ in range(5)]:
            if not 1024 * PhiH >= 2 * B * sb - B * B - sb2 - so2:
                viol["N3L"] += 1
        up = Bbar * sb - sb2 + (Mtot - 1) * so2
        if not 1024 * PhiH <= up:
            viol["N3U"] += 1
        slack = (up - 1024 * PhiH) / max(up, F(1))
        worst["N3U"] = slack if worst["N3U"] is None or slack < worst["N3U"] else worst["N3U"]
        if 1024 * PhiH > (6 * m) * sb - sb2 + (Mtot - 1) * so2:
            viol["NEGCTL"] += 1
        for t0 in t0s:
            He = (Hval(mem, t0) + Hval(memr, t0)) / 2
            d = He - F(1, 256) * t0 ** 8
            if d < 0:
                viol["POS"] += 1
            worst["POS"] = d if worst["POS"] is None or d < worst["POS"] else worst["POS"]
    print(f"  H-side cut tests (lin = {lin}, {ntr} member sets, {len(Bs)} certificate B's + 5 random, {len(t0s)} certificate t0's):")
    print(f"    violations: {viol}  (all must be 0 except NEGCTL > 0)")
    print(f"    min relative N3U slack {float(worst['N3U']):.3e}; min He(t0) - 2^-8 t0^8 = {float(worst['POS']) if worst['POS'] is not None else None}")


if __name__ == "__main__":
    for lin, fn, Sb in [(0, "xn_cert9_l0.gp", 7.3105389520668723), (1, "xn_cert9_l1.gp", 6.2989094049298407)]:
        kinds = coverage(fn, lin, Sb)
        cut_tests(lin, kinds, int(sys.argv[1]) if len(sys.argv) > 1 else 300)
