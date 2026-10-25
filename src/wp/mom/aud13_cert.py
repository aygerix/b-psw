#!/usr/bin/env python3
"""aud13: INDEPENDENT exact re-check of the xn k = 9 Farkas certificates (candidate Theorem 13).

Rows are rebuilt from aud13_rows.json (aud13_mom.py: independent e_j/p_j derivation, validated by exact polynomial
division) and from the constraint definitions re-derived in the audit:
  equalities  mass: nu(D) = m;  intv: int v = S;  E0 = Ve(-1) - 2^-8 = 0;  E1 = Ve'(-1) + 2^-5 = 0
  CS   P - (S/2) M12 <= 0                     (Cauchy-Schwarz, P = M11^2/2 >= 0)
  PHI  Phi >= 0,  Phi = Ve''(-1) - 56 2^-8    (family: Phi = (e2(beta) + e2(gamma))/2^10)
  N3L  Phi - 2^-10 [2 Bh Sb - Bh^2 - sum b^2 - sum om^2] >= 0     (b = 2 + 1/(2v), om = w/v; linear member b = 2, om = 0)
  N3U  2^-10 [Bbar Sb - sum b^2 + (Mtot - 1) sum om^2] - Phi >= 0  (Bbar = 18 m + 2 lin)
  POS  Ve(t0) - 2^-8 t0^8 >= 0
Positivity of Psi is checked by EXACT rational interval subdivision in the (w, v) plane (w = u - 1/2 >= 0 by evenness),
region {1/32 <= v <= vmax, (w + 1/2)^2 <= 8 v}, lower bounds from the exact Taylor expansion at each box centre.
Usage: aud13_cert.py certfile lin [first last]
"""
import json
import sys
import time
from fractions import Fraction as Fr
from math import isqrt

try:
    from flint import fmpq as Q

    def toQ(x):
        x = Fr(x)
        return Q(x.numerator, x.denominator)
except ImportError:  # pragma: no cover
    Q = Fr

    def toQ(x):
        return Fr(x)


# ---------------- gp vector parser ----------------
def parse_gp(sx):
    pos = 0

    def skip():
        nonlocal pos
        while pos < len(sx) and sx[pos] in " \t\r\n":
            pos += 1

    def val():
        nonlocal pos
        skip()
        if sx[pos] == "[":
            pos += 1
            out = []
            skip()
            if sx[pos] == "]":
                pos += 1
                return out
            while True:
                out.append(val())
                skip()
                if sx[pos] == ",":
                    pos += 1
                    continue
                if sx[pos] == "]":
                    pos += 1
                    return out
                raise ValueError("parse error at " + str(pos))
        if sx[pos] == '"':
            e = sx.index('"', pos + 1)
            st = sx[pos + 1:e]
            pos = e + 1
            return st
        st = pos
        while pos < len(sx) and (sx[pos].isdigit() or sx[pos] in "-/+"):
            pos += 1
        return Fr(sx[st:pos])

    v = val()
    skip()
    if pos != len(sx):
        raise ValueError("trailing text")
    return v


# ---------------- rows ----------------
ROWS = json.load(open("aud13_rows.json"))


def form(lin, name):
    d = ROWS[str(lin)][name]
    col = {}
    for k, (c0, c1) in d["col"].items():
        a, j = map(int, k.split(","))
        col[(a, j)] = [Fr(c0), Fr(c1)]
    return col, Fr(d["pc"]), [Fr(x) for x in d["c"]]


def lin_comb(terms):
    """terms: list of (factor, (col, pc, c)) -> (col, pc, c)"""
    col, pc, c = {}, Fr(0), [Fr(0)] * 3
    for f, (cl, p, cc) in terms:
        for k, (x0, x1) in cl.items():
            col.setdefault(k, [Fr(0), Fr(0)])
            col[k][0] += f * x0
            col[k][1] += f * x1
        pc += f * p
        for i in range(3):
            c[i] += f * cc[i]
    return col, pc, c


def eq_row(X):            # X = 0      ->  col.nu + pc P = -c(S)
    col, pc, c = X
    return col, pc, [-x for x in c]


def ge_row(X):            # X >= 0     ->  -col.nu - pc P <= c(S)
    col, pc, c = X
    return {k: [-a, -b] for k, (a, b) in col.items()}, -pc, list(c)


def le_row(X):            # X <= 0     ->  col.nu + pc P <= -c(S)
    return eq_row(X)


def n3_bracket_lower(Bh, lin):
    """2 Bh Sb - Bh^2 - sum b^2 - sum om^2 as a form: per-member + constant"""
    col = {(0, 0): [4 * Bh - 4, Fr(0)], (-1, 0): [Bh - 2, Fr(0)], (-2, 0): [Fr(-1, 4), Fr(0)], (-2, 1): [Fr(-1), Fr(0)]}
    return col, Fr(0), [lin * (4 * Bh - 4) - Bh * Bh, Fr(0), Fr(0)]


def n3_bracket_upper(lin):
    m = 9 - lin
    Bbar, Mtot = 18 * m + 2 * lin, m + lin
    col = {(0, 0): [Fr(2 * Bbar - 4), Fr(0)], (-1, 0): [Fr(Bbar, 2) - 2, Fr(0)], (-2, 0): [Fr(-1, 4), Fr(0)],
           (-2, 1): [Fr(Mtot - 1), Fr(0)]}
    return col, Fr(0), [Fr(lin * (2 * Bbar - 4)), Fr(0), Fr(0)]


def pos_form(t0, lin):
    terms = [(t0 ** (8 - n), form(lin, f"T{n}")) for n in range(9)]
    X = lin_comb(terms)
    X[2][0] -= t0 ** 8 / 256
    return X


def eq_rows(lin):
    m = 9 - lin
    mass = ({(0, 0): [Fr(1), Fr(0)]}, Fr(0), [Fr(m), Fr(0), Fr(0)])
    intv = ({(1, 0): [Fr(1), Fr(0)]}, Fr(0), [Fr(0), Fr(1), Fr(0)])
    return [mass, intv, eq_row(form(lin, "E0")), eq_row(form(lin, "E1"))]


def ineq_row(lin, kind, prm):
    q = Fr(1, 1024)
    PHI = form(lin, "PHI")
    if kind == "CS":
        return le_row(({(1, 1): [Fr(0), Fr(-1, 2)]}, Fr(1), [Fr(0)] * 3))
    if kind == "PHI":
        return ge_row(PHI)
    if kind == "N3L":
        return ge_row(lin_comb([(Fr(1), PHI), (-q, n3_bracket_lower(Fr(prm), lin))]))
    if kind == "N3U":
        return ge_row(lin_comb([(q, n3_bracket_upper(lin)), (Fr(-1), PHI)]))
    if kind == "POS":
        return ge_row(pos_form(Fr(prm), lin))
    raise ValueError(kind)


# ---------------- exact positivity by subdivision in (w, v) ----------------
def taylor_shift(c, x0):
    """coefficients of p(x0 + x) given p's coefficients c (low -> high)"""
    c = list(c)
    n = len(c)
    for i in range(n - 1):
        for j in range(n - 2, i - 1, -1):
            c[j] = c[j] + x0 * c[j + 1]
    return c


def box_coeffs(Qwv, wc, vc, hw, hv):
    """Qwv[i][j]: coefficient of w^i v^j.  Returns C[p][q]: coefficient of x^p y^q in Q(wc + hw x, vc + hv y)"""
    ni, nj = len(Qwv), len(Qwv[0])
    A = []                                   # A[j][p]
    for j in range(nj):
        col = taylor_shift([Qwv[i][j] for i in range(ni)], wc)
        f = Q(1)
        for p in range(ni):
            col[p] = col[p] * f
            f = f * hw
        A.append(col)
    C = [[None] * nj for _ in range(ni)]
    for p in range(ni):
        row = taylor_shift([A[j][p] for j in range(nj)], vc)
        f = Q(1)
        for qq in range(nj):
            C[p][qq] = row[qq] * f
            f = f * hv
    return C


def lower_bound(C):
    lb = C[0][0]
    dw = Q(0)
    dv = Q(0)
    for p in range(len(C)):
        for q in range(len(C[0])):
            if p == 0 and q == 0:
                continue
            x = C[p][q]
            if x == 0:
                continue
            if p % 2 == 0 and q % 2 == 0:
                if x < 0:
                    lb += x
                    if p:
                        dw += -x
                    else:
                        dv += -x
            else:
                ax = x if x > 0 else -x
                lb -= ax
                if p:
                    dw += ax
                if q:
                    dv += ax
    return lb, dw, dv


def sqrt_up(x):
    """rational upper bound of sqrt(x) for rational x >= 0"""
    x = Fr(x)
    p, q = x.numerator, x.denominator
    # sqrt(p/q) = sqrt(p q)/q
    r = isqrt(p * q)
    if r * r < p * q:
        r += 1
    return Fr(r, q)


def check_positive(Qwv, Rb, vmax, maxboxes=3_000_000):
    """exact: Q(w, v) > 0 on {0 <= w, (w + 1/2)^2/8 <= v <= vmax} (requires w <= Rb there)."""
    Qq = [[toQ(x) for x in row] for row in Qwv]
    stack = [(Fr(0), Fr(Rb), Fr(1, 32), Fr(vmax), 0)]
    nbox = 0
    minlb = None
    maxdepth = 0
    while stack:
        w0, w1, v0, v1, d = stack.pop()
        nbox += 1
        if nbox > maxboxes:
            return False, nbox, maxdepth, "box limit"
        vcurve0 = (w0 + Fr(1, 2)) ** 2 / 8
        if v1 < vcurve0:
            continue                                   # box outside the region
        v0e = max(v0, vcurve0)                         # shrink to the bounding box of box /\ region
        w1e = min(w1, sqrt_up(8 * v1) - Fr(1, 2))
        if w1e < w0:
            w1e = w0
        wc, vc = (w0 + w1e) / 2, (v0e + v1) / 2
        hw, hv = (w1e - w0) / 2, (v1 - v0e) / 2
        C = box_coeffs(Qq, toQ(wc), toQ(vc), toQ(hw), toQ(hv))
        lb, dw, dv = lower_bound(C)
        if lb > 0:
            continue
        if C[0][0] <= 0:
            # the centre itself: is it in the region?  (w + 1/2)^2 <= 8 v
            if (wc + Fr(1, 2)) ** 2 <= 8 * vc:
                return False, nbox, maxdepth, f"Q <= 0 at region point (w, v) = ({float(wc)}, {float(vc)}): {float(Fr(str(C[0][0])))}"
        if d > 80:
            return False, nbox, maxdepth, f"depth limit near (w, v) = ({float(wc)}, {float(vc)})"
        maxdepth = max(maxdepth, d + 1)
        if dw >= dv:
            wm = (w0 + w1) / 2
            stack.append((w0, wm, v0, v1, d + 1))
            stack.append((wm, w1, v0, v1, d + 1))
        else:
            vm = (v0 + v1) / 2
            stack.append((w0, w1, v0, vm, d + 1))
            stack.append((w0, w1, vm, v1, d + 1))
    return True, nbox, maxdepth, "ok"


def psi_poly(col, S):
    """Q(w, v) = v^2 Psi_S(v, W = w^2): Qwv[i][j] coefficient of w^i v^j"""
    Qwv = [[Fr(0)] * 6 for _ in range(5)]
    for (a, j), (c0, c1) in col.items():
        cf = c0 + c1 * S
        if cf == 0:
            continue
        e = a + 2
        if not (0 <= e <= 5 and 0 <= 2 * j <= 4):
            raise ValueError(f"basis term v^{a} W^{j} outside the expected range")
        Qwv[2 * j][e] += cf
    return Qwv


def SofR(R, m):
    return (R + Fr(1, 2)) ** 2 / 8 + Fr(m - 1, 32)


def check_cert(c):
    lin, Ra, Rb, y, zl = c
    lin = int(lin)
    m = 9 - lin
    Sa, Sb = SofR(Ra, m), SofR(Rb, m)
    E = eq_rows(lin)
    if len(y) != 4:
        return False, "y length"
    terms = [(y[i], E[i]) for i in range(4)]
    for kind, prm, z in zl:
        if z < 0:
            return False, "negative z"
        terms.append((z, ineq_row(lin, kind, prm)))
    col, pc, kap = lin_comb(terms)
    if pc < 0:
        return False, f"P column negative ({float(pc)})"
    cand = [Sa, Sb]
    if kap[2] != 0:
        sv = -kap[1] / (2 * kap[2])
        if Sa < sv < Sb:
            cand.append(sv)
    kmax = max(kap[0] + kap[1] * x + kap[2] * x * x for x in cand)
    if kmax >= 0:
        return False, f"kappa not negative (max {float(kmax)})"
    vmax = Sb - Fr(m - 1, 32)
    assert vmax == (Rb + Fr(1, 2)) ** 2 / 8
    info = []
    for Sx in (Sa, Sb):
        Qwv = psi_poly(col, Sx)
        ok, nb, md, msg = check_positive(Qwv, Rb, vmax)
        info.append(f"S={float(Sx):.5f}: boxes {nb} depth {md}")
        if not ok:
            return False, f"Psi positivity FAILED at S = {float(Sx)}: {msg} ({nb} boxes)"
    used = sorted(set(k for k, _, _ in zl))
    return True, f"pc {float(pc):.3g}, kappa max {float(kmax):.4g}, ineqs {used} x{len(zl)}; " + "; ".join(info)


def main():
    fname, lin = sys.argv[1], int(sys.argv[2])
    L = [parse_gp(line.strip()) for line in open(fname) if line.strip()]
    first = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    last = int(sys.argv[4]) if len(sys.argv) > 4 else len(L)
    bad = 0
    for i in range(first - 1, last):
        c = L[i]
        if int(c[0]) != lin:
            print(f"cert {i + 1}: lin mismatch")
            bad += 1
            continue
        t = time.time()
        ok, msg = check_cert(c)
        bad += 0 if ok else 1
        m = 9 - lin
        print(f"cert {i + 1:3d} lin={lin} R=[{c[1]}, {c[2]}] S=[{float(SofR(c[1], m)):.6f}, {float(SofR(c[2], m)):.6f}] "
              f"{'PASS' if ok else 'FAIL'} ({time.time() - t:.1f}s): {msg}", flush=True)
    print(f"{fname}: certificates {first}..{last}: {bad} failed", flush=True)


if __name__ == "__main__":
    main()
