#!/usr/bin/env python3
"""aud13 controls for aud13_cert.py.
 (1) sharpness: for every certificate, estimate mu = min Psi over the region (float grid, both S endpoints), then shift the
     mass multiplier y_1 by -1.01 mu (Psi -> Psi - 1.01 mu, kappa more negative): must FAIL; by -0.5 mu: should PASS.
 (2) corrupted rows: PHI coefficient of M2_0 + 1; POS with t0 + 1/5; E0 constant * (1 + 1e-6): count failures.
 (3) wrong case: lin = 0 certificates checked with lin = 1 rows.
 (4) coverage: union of certified R-intervals vs the S-bound (exact).
"""
import copy
import sys
from fractions import Fraction as Fr

import numpy as np
import sympy as sp

import aud13_cert as AC


def psi_min_grid(c, n_r=400, n_t=41):
    lin, Ra, Rb, y, zl = c
    lin = int(lin)
    m = 9 - lin
    E = AC.eq_rows(lin)
    terms = [(y[i], E[i]) for i in range(4)] + [(z, AC.ineq_row(lin, k, p)) for k, p, z in zl]
    col, pc, kap = AC.lin_comb(terms)
    mus = []
    rs = np.linspace(0, float(Rb), n_r)
    for Sx in (AC.SofR(Ra, m), AC.SofR(Rb, m)):
        vals = None
        for th in np.linspace(0, 1, n_t):
            v = (rs + 0.5) ** 2 / 8
            W = th * rs ** 2
            tot = np.zeros_like(rs)
            for (a, j), (c0, c1) in col.items():
                tot += float(c0 + c1 * Sx) * v ** a * W ** j
            mn = tot.min()
            vals = mn if vals is None else min(vals, mn)
        mus.append(vals)
    return min(mus)


def main():
    certs = {0: [AC.parse_gp(l.strip()) for l in open("xn_cert9_l0.gp") if l.strip()],
             1: [AC.parse_gp(l.strip()) for l in open("xn_cert9_l1.gp") if l.strip()]}
    # (1) sharpness
    nf_hi = nf_lo = 0
    for lin in (0, 1):
        for i, c in enumerate(certs[lin]):
            mu = psi_min_grid(c)
            res = []
            for fac in (1.01, 0.5):
                cc = copy.deepcopy(c)
                cc[3][0] = cc[3][0] - Fr(fac * mu).limit_denominator(10 ** 12)
                ok, msg = AC.check_cert(cc)
                res.append(ok)
            nf_hi += 0 if res[0] else 1
            nf_lo += 0 if res[1] else 1
            print(f"ctl1 lin={lin} cert {i + 1:2d}: grid min Psi = {mu:.6g}; shift 1.01 mu -> {'PASS (BAD)' if res[0] else 'FAIL (good)'}; "
                  f"shift 0.5 mu -> {'PASS' if res[1] else 'FAIL'}", flush=True)
    print(f"ctl1 summary: shifted by 1.01 mu rejected {nf_hi}/109 (expect 109); shifted by 0.5 mu rejected {109 - (109 - nf_lo)}/109 (expect ~0)")
    # (2) corrupted rows
    orig = copy.deepcopy(AC.ROWS)
    for label, mut in (("PHI M2_0 + 1", lambda R: R[str(l)]["PHI"]["col"].__setitem__("2,0", [str(Fr(R[str(l)]["PHI"]["col"]["2,0"][0]) + 1), R[str(l)]["PHI"]["col"]["2,0"][1]])),
                       ("E0 const x (1 + 1e-6)", lambda R: R[str(l)]["E0"]["c"].__setitem__(0, str(Fr(R[str(l)]["E0"]["c"][0]) * (1 + Fr(1, 10 ** 6))))),
                       ("T5 M1_2 coefficient + 1/100", lambda R: R[str(l)]["T5"]["col"].__setitem__("1,1", [str(Fr(R[str(l)]["T5"]["col"]["1,1"][0]) + Fr(1, 100)), R[str(l)]["T5"]["col"]["1,1"][1]]))):
        tot = fails = 0
        for l in (0, 1):
            AC.ROWS = copy.deepcopy(orig)
            mut(AC.ROWS)
            for c in certs[l]:
                tot += 1
                ok, _ = AC.check_cert(c)
                fails += 0 if ok else 1
        print(f"ctl2 corrupted rows [{label}]: {fails}/{tot} certificates rejected", flush=True)
    AC.ROWS = orig
    # POS t0 shifted
    fails = tot = 0
    for l in (0, 1):
        for c in certs[l]:
            if not any(k == "POS" for k, _, _ in c[4]):
                continue
            cc = copy.deepcopy(c)
            cc[4] = [[k, (p + Fr(1, 5)) if k == "POS" else p, z] for k, p, z in cc[4]]
            tot += 1
            ok, _ = AC.check_cert(cc)
            fails += 0 if ok else 1
    print(f"ctl2 POS t0 + 1/5: {fails}/{tot} POS-using certificates rejected", flush=True)
    # (3) wrong case
    fails = 0
    for c in certs[0]:
        cc = copy.deepcopy(c)
        cc[0] = Fr(1)
        ok, _ = AC.check_cert(cc)
        fails += 0 if ok else 1
    print(f"ctl3 lin=0 certificates with lin=1 rows: {fails}/57 rejected", flush=True)
    # (4) coverage
    for lin, lam, rhs_c, phimax in ((0, Fr(59, 110), (Fr(39391023, 1100000), Fr(-788, 125), Fr(-2, 55)), Fr(729, 32)),
                                    (1, Fr(27, 50), (Fr(1294023, 50000), Fr(-678, 125), Fr(-1, 25)), Fr(585, 32))):
        m = 9 - lin
        iv = sorted((c[1], c[2]) for c in certs[lin])
        reach = Fr(0)
        for a, b in iv:
            if a > reach:
                break
            reach = max(reach, b)
        Srch = AC.SofR(reach, m)
        S = sp.Symbol("S")
        poly = sp.Rational(rhs_c[0] + lam * phimax) + sp.Rational(rhs_c[1]) * S + sp.Rational(rhs_c[2]) * S ** 2
        root = max(r for r in sp.solve(poly, S) if r.is_real)
        val = poly.subs(S, sp.Rational(Srch))
        print(f"ctl4 lin={lin}: certified R-intervals are contiguous from 0 to {reach}; S covered (m/32, {float(Srch):.6f}]; "
              f"S-bound root {sp.N(root, 10)}; bound polynomial at the cover end = {sp.N(val, 6)} ({'< 0: COVERED' if val < 0 else 'NOT covered'})")


if __name__ == "__main__":
    main()
