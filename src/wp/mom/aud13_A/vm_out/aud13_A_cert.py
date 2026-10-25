#!/usr/bin/env python3
"""aud13_A_cert.py -- audit of candidate Theorem 13: re-check the xn Farkas certificates with MY rows and a different
positivity test.  New code; reads xn_cert9_l{0,1}.gp only as data (lin, Ra, Rb, y, [kind, param, z]).

Rows (semantics: column . nu + pc P  (= or <=)  rhs(S); nu = counting measure of the cubic members on the closed
domain D_{Sb}; P = M11^2/2):
  equalities  mass (col 1, rhs m), S (col v, rhs S), Z0e = Ve(-1) - 2^-8 = 0, Z1e = Ve'(-1) + 8 2^-8 = 0
  CS   P - (S/2) sum v W <= 0                      PHI  -Phi_V <= 0
  POS  -(Ve(t0) - 2^-8 t0^8) <= 0                  N3L/N3U implemented for completeness (b = 2 + 1/(2v), om^2 = W/v^2)
Moment forms come from aud13_A_rows.txt (my derivation, aud13_A_mom.gp).
Positivity: Q = v^2 Psi_S(W, v) is a polynomial in (W, v); the closed domain is v in [1/32, vcap], 0 <= W <= Wmax(v),
Wmax(v) = (sqrt(8v) - 1/2)^2.  Boxes in (W, v) (NOT the (r, theta) chart of the line), exact rational Bernstein bounds,
boxes trimmed with rigorous rational bounds of the curved boundary; a box is accepted iff all its Bernstein
coefficients are > 0, discarded iff it lies outside the domain (exact test), otherwise split.
Usage: aud13_A_cert.py lin first last [maxboxes]
"""
import re
import sys
import time
from fractions import Fraction as F
from math import comb, isqrt

LINMONO = [(1, 1), (2, 0), (1, 2), (2, 1), (3, 0)]     # M1_2, M2_0, M1_4, M2_2, M3_0 per member: v W, v^2, v W^2, v^2 W, v^3


def load_forms(fn="aud13_A_rows.txt"):
    out = {}
    for line in open(fn):
        name, lin, rest = line.split(" ", 2)
        parts = rest.split("|")
        c = [F(x) for x in parts[1].split()]
        d = [F(x) for x in parts[2].split()]
        q = F(parts[3].strip())
        col = {}
        for i, mono in enumerate(LINMONO):
            if d[2 * i] or d[2 * i + 1]:
                col[mono] = [d[2 * i], d[2 * i + 1]]
        out[(name, int(lin))] = dict(col=col, pc=2 * q, c=c)      # X = col.nu + pc P + c(S)
    return out


def addcol(dst, src, f):
    for k, (a, b) in src.items():
        x = dst.setdefault(k, [F(0), F(0)])
        x[0] += f * a
        x[1] += f * b


def row(col, pc, rhs):
    return dict(col=col, pc=F(pc), rhs=[F(x) for x in rhs])


def build_rows(FM, lin, kinds):
    m = 9 - lin
    eq = [row({(0, 0): [F(1), F(0)]}, 0, [m, 0, 0]), row({(1, 0): [F(1), F(0)]}, 0, [0, 1, 0])]
    for nm in ("Z0", "Z1"):
        X = FM[(nm, lin)]
        eq.append(row({k: list(v) for k, v in X["col"].items()}, X["pc"], [-x for x in X["c"]]))
    ph = FM[("PHI", lin)]
    ineq = []
    for (kd, prm, z) in kinds:
        if kd == "CS":
            r = row({(1, 1): [F(0), F(-1, 2)]}, 1, [0, 0, 0])
        elif kd == "PHI":
            col = {}
            addcol(col, ph["col"], -1)
            r = row(col, -ph["pc"], ph["c"])
        elif kd == "POS":
            t0 = F(prm)
            col, pc, rhs = {}, F(0), [(1 - F(1, 256)) * t0 ** 8, F(0), F(0)]
            for n in range(1, 9):
                X = FM[(f"T{n}", lin)]
                p = t0 ** (8 - n)
                addcol(col, X["col"], -p)
                pc -= p * X["pc"]
                for i in range(3):
                    rhs[i] += p * X["c"][i]
            r = row(col, pc, rhs)
        elif kd in ("N3L", "N3U"):
            Mtot, Bbar = m + lin, 18 * m + 2 * lin
            bcol = {(0, 0): [F(2), F(0)], (-1, 0): [F(1, 2), F(0)]}
            b2col = {(0, 0): [F(4), F(0)], (-1, 0): [F(2), F(0)], (-2, 0): [F(1, 4), F(0)]}
            o2col = {(-2, 1): [F(1), F(0)]}
            col = {}
            if kd == "N3L":
                B = F(prm)
                addcol(col, ph["col"], -1)
                addcol(col, bcol, 2 * B / 1024); addcol(col, b2col, F(-1, 1024)); addcol(col, o2col, F(-1, 1024))
                linc = (2 * B * 2 - 4) * lin - B * B          # lin member's 2B b - b^2 - om^2, and -B^2
                r = row(col, -ph["pc"], [ph["c"][0] - linc / 1024, ph["c"][1], ph["c"][2]])
            else:
                addcol(col, ph["col"], 1)
                addcol(col, bcol, F(-Bbar, 1024)); addcol(col, b2col, F(1, 1024)); addcol(col, o2col, F(-(Mtot - 1), 1024))
                linc = (Bbar * 2 - 4) * lin
                r = row(col, ph["pc"], [-ph["c"][0] + linc / 1024, -ph["c"][1], -ph["c"][2]])
        else:
            raise ValueError(kd)
        ineq.append((r, z))
    return eq, ineq


# ---------------- positivity on the curved domain ----------------
SQD = 2 ** 64


def sqrt_lo(x):
    x = F(x)
    return F(isqrt(x.numerator * x.denominator * SQD * SQD), x.denominator * SQD)


def sqrt_hi(x):
    return sqrt_lo(x) + F(1, F(x).denominator * SQD)


DY = 2 ** 48


def floor_dy(x):
    return F((x.numerator * DY) // x.denominator, DY)


def ceil_dy(x):
    return F(-((-x.numerator * DY) // x.denominator), DY)


def outside(W0, v1):
    """exact: is W0 > Wmax(v1), i.e. (sqrt(W0) + 1/2)^2 > 8 v1 ?"""
    rhs = 8 * v1 - W0 - F(1, 4)          # need sqrt(W0) > rhs
    if rhs < 0:
        return True
    return W0 > rhs * rhs


def shift2(A, n1, n2, W0, hW, v0, hv):
    """A[i][j] coefficients of W^i v^j -> coefficients of x^i y^j with W = W0 + hW x, v = v0 + hv y"""
    B = [[F(0)] * (n2 + 1) for _ in range(n1 + 1)]
    pw = [[F(comb(i, a)) * W0 ** (i - a) * hW ** a for a in range(i + 1)] for i in range(n1 + 1)]
    pv = [[F(comb(j, b)) * v0 ** (j - b) * hv ** b for b in range(j + 1)] for j in range(n2 + 1)]
    for i in range(n1 + 1):
        for j in range(n2 + 1):
            a_ = A[i][j]
            if a_ == 0:
                continue
            for a in range(i + 1):
                ca = a_ * pw[i][a]
                for b in range(j + 1):
                    B[a][b] += ca * pv[j][b]
    return B


def bern_min(Bc, n1, n2):
    mn, arg = None, None
    for k in range(n1 + 1):
        for l in range(n2 + 1):
            s = F(0)
            for i in range(k + 1):
                ci = F(comb(k, i), comb(n1, i))
                for j in range(l + 1):
                    if Bc[i][j]:
                        s += ci * F(comb(l, j), comb(n2, j)) * Bc[i][j]
            if mn is None or s < mn:
                mn, arg = s, (k, l)
    return mn, arg


def evalQ(A, W, v):
    return sum(A[i][j] * W ** i * v ** j for i in range(len(A)) for j in range(len(A[0])))


def check_pos(A, vcap, maxboxes):
    """A: coefficients of Q(W, v) (degree n1 in W, n2 in v).  Returns (ok, nboxes, info)."""
    n1, n2 = len(A) - 1, len(A[0]) - 1
    v_lo = F(1, 32)
    W_hi = ceil_dy((sqrt_hi(8 * vcap) - F(1, 2)) ** 2)
    stack = [(F(0), W_hi, v_lo, vcap, 0)]
    nb, maxd = 0, 0
    while stack:
        W0, W1, v0, v1, d = stack.pop()
        if outside(W0, v1):
            continue
        # trim: W <= Wmax(v1) <= W1', v >= vmin(W0) >= v0'
        W1 = min(W1, ceil_dy((sqrt_hi(8 * v1) - F(1, 2)) ** 2))
        v0 = max(v0, floor_dy((sqrt_lo(W0) + F(1, 2)) ** 2 / 8))
        if W1 < W0 or v1 < v0:
            continue
        nb += 1
        maxd = max(maxd, d)
        if nb > maxboxes:
            return False, nb, f"box budget exceeded (depth {maxd})"
        hW, hv = W1 - W0, v1 - v0
        Bc = shift2(A, n1, n2, W0, hW, v0, hv)
        mn, arg = bern_min(Bc, n1, n2)
        if mn > 0:
            continue
        # genuine negative at a domain corner?
        for (Wc, vc) in ((W0, v0), (W0, v1), (W1, v0), (W1, v1)):
            if not outside(Wc, vc) and Wc <= (sqrt_lo(8 * vc) - F(1, 2)) ** 2 and evalQ(A, Wc, vc) <= 0:
                return False, nb, f"Q <= 0 at domain point W = {float(Wc):.6g}, v = {float(vc):.6g}: {float(evalQ(A, Wc, vc)):.3g}"
        if hW == 0 and hv == 0:
            return False, nb, "degenerate box"
        # split the relatively larger side (scale W by its full range, v by its range)
        if hW * (vcap - F(1, 32)) >= hv * W_hi:
            Wm = (W0 + W1) / 2
            stack.append((W0, Wm, v0, v1, d + 1)); stack.append((Wm, W1, v0, v1, d + 1))
        else:
            vm = (v0 + v1) / 2
            stack.append((W0, W1, v0, vm, d + 1)); stack.append((W0, W1, vm, v1, d + 1))
    return True, nb, f"depth {maxd}"


def Qcoef(col, S):
    """v^2 * Psi_S as A[i][j] (W^i v^j)"""
    terms = {}
    for (p, q), (a, b) in col.items():
        c = a + b * S
        if c:
            terms[(q, p + 2)] = terms.get((q, p + 2), F(0)) + c
    assert all(j >= 0 for (_, j) in terms)
    n1 = max([i for (i, _) in terms] + [0])
    n2 = max([j for (_, j) in terms] + [0])
    A = [[F(0)] * (n2 + 1) for _ in range(n1 + 1)]
    for (i, j), c in terms.items():
        A[i][j] = c
    return A


def SofR(R, m):
    return (R + F(1, 2)) ** 2 / 8 + F(m - 1, 32)


def check_cert(FM, c, maxboxes):
    lin, Ra, Rb, y, zl = c
    m = 9 - lin
    Sa, Sb = SofR(Ra, m), SofR(Rb, m)
    eq, ineq = build_rows(FM, lin, zl)
    col, pc, kap = {}, F(0), [F(0)] * 3
    for yi, r in zip(y, eq):
        addcol(col, r["col"], yi); pc += yi * r["pc"]
        kap = [kap[i] + yi * r["rhs"][i] for i in range(3)]
    for (r, z) in ineq:
        if z < 0:
            return False, "negative z"
        addcol(col, r["col"], z); pc += z * r["pc"]
        kap = [kap[i] + z * r["rhs"][i] for i in range(3)]
    if pc < 0:
        return False, f"P column {pc} < 0"
    cand = [Sa, Sb]
    if kap[2] != 0:
        sv = -kap[1] / (2 * kap[2])
        if Sa < sv < Sb:
            cand.append(sv)
    kmax = max(kap[0] + kap[1] * s + kap[2] * s * s for s in cand)
    if kmax >= 0:
        return False, f"kappa max {float(kmax)} >= 0"
    vcap = Sb - F(m - 1, 32)
    info = []
    for S in (Sa, Sb):
        A = Qcoef(col, S)
        ok, nb, msg = check_pos(A, vcap, maxboxes)
        info.append(f"S={float(S):.4f}: {nb} boxes {msg}")
        if not ok:
            return False, "positivity: " + "; ".join(info)
    return True, f"pc {float(pc):.3g}, kappa max {float(kmax):.4g}; " + "; ".join(info)


def load_certs(fn):
    out = []
    for line in open(fn):
        line = line.strip()
        if line:
            out.append(eval(re.sub(r"(-?\d+)/(\d+)", r"F(\1,\2)", line), {"__builtins__": {}}, {"F": F}))
    return out


if __name__ == "__main__":
    lin = int(sys.argv[1])
    first, last = int(sys.argv[2]), int(sys.argv[3])
    maxboxes = int(sys.argv[4]) if len(sys.argv) > 4 else 200000
    mode = sys.argv[5] if len(sys.argv) > 5 else "none"     # negative controls: dropPHI, dropPOS, dropCS, scaleZ0
    FM = load_forms()
    C = load_certs(f"xn_cert9_l{lin}.gp")
    if mode.startswith("drop"):
        kd = mode[4:]
        C = [[c[0], c[1], c[2], c[3], [x for x in c[4] if x[0] != kd]] for c in C]
    elif mode == "scaleZ0":
        C = [[c[0], c[1], c[2], [c[3][0], c[3][1], c[3][2] * F(101, 100), c[3][3]], c[4]] for c in C]
    nfail = 0
    for idx in range(first, min(last, len(C)) + 1):
        c = C[idx - 1]
        t0 = time.time()
        ok, msg = check_cert(FM, c, maxboxes)
        nfail += (not ok)
        kinds = sorted(set(k + ("" if k != "POS" else str(p)) for k, p, z in c[4]))
        print(f"cert {idx:3d} lin={lin} R=[{c[1]}, {c[2]}] {'PASS' if ok else 'FAIL'} ({time.time() - t0:.1f}s) rows {kinds}: {msg}", flush=True)
    print(f"lin={lin} certs {first}..{min(last, len(C))}: {nfail} failed")
