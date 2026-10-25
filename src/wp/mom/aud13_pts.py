#!/usr/bin/env python3
"""aud13 (a): evaluate EVERY xn relaxation constraint at every m9 point, k = 9.

Points (reconstructed from the m9 data files, then re-solved here with the DIRECT equations):
  X1..X4  the exact symmetric level-1 points of m9_symexact.out (pair choices [1,2], [1,3], [1,4], [2,3]), 60 digits;
  C90     certified level-1 point, 9 cubics (m9_kcert_9_0.out: fixed rationals + Newton on the 4 unknowns);
  C81     certified level-1 point, 8 cubics + 2t+1 (m9_kcert_9_1.out);
  L90     certified level-2 point, 9 cubics (m9_l2cert_9.gp start, fixed = round(x 10^8)/10^8, Newton on 6 unknowns);
  L81     certified level-2 point, 8 cubics + 2t+1.
Direct quantities: V = (F - 1) div (t prod g) by polynomial division in mpmath (no moment forms), V~ for the
reflected family, beta_i = 2 + (1-u)/v, gamma_i = 2 + u/v (= 2 for 2t + 1).
Constraints: E0, E1 (and the odd parts), CS, PHI, N3L (all Bh, i.e. Phi >= 2^-10 (Sb^2 - sum b^2 - sum om^2)), N3U,
POS (t0 list of xn plus a fine grid), domain (v > 1/32, |w| < sqrt(8v) - 1/2), S-identity at t0 = -7/5, S-bound.
For every certificate whose S-interval contains S: residual of every row it uses (DIRECT values) and the Farkas identity
  sum_members Psi + pc P - kappa(S) = sum_i y_i res_i + sum_j z_j res_j   (left side from the aud13 moment-form rows).
"""
import re
from fractions import Fraction as Fr

import mpmath as mp

import aud13_cert as AC

mp.mp.dps = 70
K = 9


# ---------------- direct polynomial computations ----------------
def pmul(a, b):
    o = [mp.mpf(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            o[i + j] += x * y
    return o


def pquo(num, den):
    num = list(num)
    q = [mp.mpf(0)] * (len(num) - len(den) + 1)
    for sh in range(len(num) - len(den), -1, -1):
        f = num[sh + len(den) - 1] / den[-1]
        q[sh] = f
        for i, d in enumerate(den):
            num[sh + i] -= f * d
    return q


def Vpoly(U, Vv, lin):
    """quotient of (F - 1) by t prod g; members f = a t^3 + (a + c) t^2 + (c + 2) t + 1, a = 1/v, c = u/v"""
    F = [mp.mpf(1)]
    G = [mp.mpf(0), mp.mpf(1)]
    for u, v in zip(U, Vv):
        a, c = 1 / v, u / v
        F = pmul(F, [mp.mpf(1), c + 2, a + c, a])
        G = pmul(G, [c + 2, a + c, a])
    for _ in range(lin):
        F = pmul(F, [mp.mpf(1), mp.mpf(2)])
        G = pmul(G, [mp.mpf(2)])
    F[0] -= 1
    return pquo(F, G)


def peval(p, x, d=0):
    out = mp.mpf(0)
    for n in range(len(p) - 1, d - 1, -1):
        cf = p[n]
        for j in range(d):
            cf *= (n - j)
        out = out * x + cf
    return out


def e2(xs):
    s1 = sum(xs)
    return (s1 * s1 - sum(x * x for x in xs)) / 2


def boundary_system(U, Vv, lin, level):
    V = Vpoly(U, Vv, lin)
    Vr = Vpoly([1 - u for u in U], Vv, lin)
    r = [peval(V, -1) - mp.mpf(1) / 256, peval(V, -1, 1) + mp.mpf(8) / 256,
         peval(Vr, -1) - mp.mpf(1) / 256, peval(Vr, -1, 1) + mp.mpf(8) / 256]
    if level >= 2:
        beta = [2 + (1 - u) / v for u, v in zip(U, Vv)] + [mp.mpf(2)] * lin
        gamma = [2 + u / v for u, v in zip(U, Vv)] + [mp.mpf(2)] * lin
        r.append(peval(V, -1, 2) - mp.mpf(56) / 256 - e2(beta) / 512)
        r.append(peval(Vr, -1, 2) - mp.mpf(56) / 256 - e2(gamma) / 512)
    return r


def newton(X, unk, m, lin, level, it=60):
    X = [mp.mpf(x) for x in X]
    for _ in range(it):
        U, Vv = X[:m], X[m:]
        F = boundary_system(U, Vv, lin, level)
        nrm = max(abs(f) for f in F)
        if nrm < mp.mpf(10) ** (-60):
            break
        J = mp.matrix(len(F), len(unk))
        h = mp.mpf(10) ** (-30)
        for c, idx in enumerate(unk):
            Xp = list(X)
            Xp[idx] += h
            Xm = list(X)
            Xm[idx] -= h
            Fp = boundary_system(Xp[:m], Xp[m:], lin, level)
            Fm = boundary_system(Xm[:m], Xm[m:], lin, level)
            for r in range(len(F)):
                J[r, c] = (Fp[r] - Fm[r]) / (2 * h)
        dx = mp.lu_solve(J, mp.matrix(F))
        for c, idx in enumerate(unk):
            X[idx] -= dx[c]
    U, Vv = X[:m], X[m:]
    return X, max(abs(f) for f in boundary_system(U, Vv, lin, level))


# ---------------- points ----------------
def Q(st):
    return mp.mpf(Fr(st).numerator) / Fr(st).denominator


def sym_points():
    txt = open("m9_symexact.out").read()
    pts = []
    for blk in txt.split("==== unknown pairs ")[1:]:
        if "EXPLICIT ALGEBRAIC FEASIBLE POINT" not in blk:
            continue
        i1, i2 = map(int, re.match(r"\[(\d+), (\d+)\]", blk).groups())
        PV = [Q(x) for x in re.search(r"v \(pairs\) = \[([^\]]*)\]", blk).group(1).split(",")]
        vc = Q(re.search(r"v_c = ([0-9/]+)", blk).group(1))
        fixed = re.search(r"fixed y = \[(.*)\]\s*\n", blk).group(1)
        fy = {int(a): Q(b) for a, b in re.findall(r"\[(\d+), ([0-9/]+)\]", fixed)}
        mt = re.search(r"y_\d+ = \(([-0-9/]+) \+ \(([-0-9]+)\) sqrt\(D\)\) / \(([-0-9/]+)\),\s+y_\d+ = ([-0-9/]+) \+ \(([-0-9/]+)\) y_\d+", blk)
        A, e, B, c0, c1 = mt.groups()
        D = Q(re.search(r"\n\s+D = ([0-9/]+)", blk).group(1))
        y1 = (Q(A) + int(e) * mp.sqrt(D)) / Q(B)
        y2 = Q(c0) + Q(c1) * y1
        Y = dict(fy)
        Y[i1], Y[i2] = y1, y2
        U = [mp.mpf(1) / 2 + mp.sqrt(Y[j]) for j in range(1, 5)] + [mp.mpf(1) / 2 - mp.sqrt(Y[j]) for j in range(1, 5)] + [mp.mpf(1) / 2]
        Vv = PV + PV + [vc]
        pts.append((f"X[{i1},{i2}]", U, Vv, 0, 1))
    return pts


def kcert_point(fname, name, m, lin):
    txt = open(fname).read()
    unk = re.search(r"unknown coordinates: \[([^\]]*)\]", txt).group(1).replace(" ", "").split(",")
    fixed = [Q(x) for x in re.search(r"fixed coordinates \(rational, denominator 10\^8\): \[([^\]]*)\]", txt).group(1).split(",")]
    cen = [mp.mpf(x) for x in re.search(r"box centre of the unknowns[^:]*: \[([^\]]*)\]", txt).group(1).split(",")]
    names = [f"u{i}" for i in range(1, m + 1)] + [f"v{i}" for i in range(1, m + 1)]
    uidx = [names.index(x) for x in unk]
    X = [None] * (2 * m)
    fi = iter(fixed)
    for i in range(2 * m):
        if i not in uidx:
            X[i] = next(fi)
    for i, c in zip(uidx, cen):
        X[i] = c
    X, res = newton(X, uidx, m, lin, 1)
    corr = max(abs(X[i] - c) for i, c in zip(uidx, cen))
    print(f"  {name}: Newton residual {mp.nstr(res, 3)}, |x - printed centre| = {mp.nstr(corr, 3)} (centre printed to 30 digits)")
    return (name, X[:m], X[m:], lin, 1)


L2 = {(9, 0): ([-1.5022143532841334, -4.123908001403858, 0.5139699491914088, -6.9488109387036845, 0.4557359591171056, 1.9898717712453684,
                8.384115532862122, 0.0300074710797625, 0.2116870195358351],
               [2.1953624588934777, 4.160359143007728, 0.1638235274005706, 12.93704174320486, 0.1504054295230128, 11.332431645079225,
                13.503906596594524, 0.8344275287407846, 16.583303646561404], ["u4", "u7", "u9", "v4", "v7", "v9"]),
      (8, 1): ([7.383602524018461, -0.281220887363323, 3.7362686543302996, -6.113476957120817, -1.9894291170655451, 0.9158133556559549,
                -1.4477281105089814, -2.7740334229151604],
               [10.187126685349437, 7.09708766558695, 3.420275946834393, 9.28262974093801, 1.3912245328824524, 13.036011292752544,
                1.249721086956731, 4.659440794003997], ["u1", "u4", "u6", "v1", "v4", "v6"])}


def l2_point(m, lin):
    U0, V0, unk = L2[(m, lin)]
    names = [f"u{i}" for i in range(1, m + 1)] + [f"v{i}" for i in range(1, m + 1)]
    uidx = [names.index(x) for x in unk]
    X0 = U0 + V0
    X = []
    for i, x in enumerate(X0):
        fx = Fr(repr(x))
        X.append(mp.mpf(x) if i in uidx else Q(Fr(round(fx * 10 ** 8), 10 ** 8)))
    X, res = newton(X, uidx, m, lin, 2)
    name = f"L{m}{lin}"
    print(f"  {name}: level-2 Newton residual {mp.nstr(res, 3)}, max |x - start| = {mp.nstr(max(abs(X[i] - X0[i]) for i in uidx), 3)}")
    return (name, X[:m], X[m:], lin, 2)


# ---------------- constraint evaluation ----------------
T0LIST = [Fr(-4) + Fr(i, 5) for i in range(36)]


def evaluate(pt, certs):
    name, U, Vv, lin, level = pt
    m = len(U)
    Wd = [u - mp.mpf(1) / 2 for u in U]
    S = sum(Vv)
    M12 = sum(v * w * w for v, w in zip(Vv, Wd))
    M20 = sum(v * v for v in Vv)
    M11 = sum(v * w for v, w in zip(Vv, Wd))
    P = M11 ** 2 / 2
    V = Vpoly(U, Vv, lin)
    Vr = Vpoly([1 - u for u in U], Vv, lin)

    def Ve(t, d=0):
        return (peval(V, t, d) + peval(Vr, t, d)) / 2

    E0 = Ve(-1) - mp.mpf(1) / 256
    E1 = Ve(-1, 1) + mp.mpf(8) / 256
    O0 = (peval(V, -1) - peval(Vr, -1)) / 2
    O1 = (peval(V, -1, 1) - peval(Vr, -1, 1)) / 2
    Phi = Ve(-1, 2) - mp.mpf(56) / 256
    b = [2 + 1 / (2 * v) for v in Vv] + [mp.mpf(2)] * lin
    om = [w / v for w, v in zip(Wd, Vv)] + [mp.mpf(0)] * lin
    beta = [x - y for x, y in zip(b, om)]
    gamma = [x + y for x, y in zip(b, om)]
    N3 = (e2(beta) + e2(gamma)) / 1024
    Sb, So = sum(b), sum(om)
    sb2, so2 = sum(x * x for x in b), sum(x * x for x in om)
    N3Lbest = (Sb ** 2 - sb2 - so2) / 1024
    Bbar, Mtot = 18 * m + 2 * lin, m + lin
    N3Uval = (Bbar * Sb - sb2 + (Mtot - 1) * so2) / 1024
    pos = {t0: Ve(mp.mpf(t0.numerator) / t0.denominator) - (mp.mpf(t0.numerator) / t0.denominator) ** 8 / 256 for t0 in T0LIST}
    posmin_t0 = min(pos, key=lambda t: pos[t])
    grid = [mp.mpf(-12) + mp.mpf(i) / 100 for i in range(2401)]
    posg = [(Ve(t) - t ** 8 / 256, t) for t in grid]
    posgmin = min(posg)
    t75 = mp.mpf(-7) / 5
    Re75 = (Ve(t75) - t75 ** 8 / 256) / (t75 + 1) ** 2
    if lin == 0:
        lam, RHS, phimax = Fr(59, 110), -mp.mpf(2) / 55 * S ** 2 - mp.mpf(788) / 125 * S - mp.mpf(2) / 55 * M20 + mp.mpf(39391023) / 1100000, mp.mpf(729) / 32
    else:
        lam, RHS, phimax = Fr(27, 50), -mp.mpf(1) / 25 * S ** 2 - mp.mpf(678) / 125 * S - mp.mpf(1) / 25 * M20 + mp.mpf(1294023) / 50000, mp.mpf(585) / 32
    lamf = mp.mpf(lam.numerator) / lam.denominator
    dommarg = min(min(8 * v - (abs(w) + mp.mpf(1) / 2) ** 2 for v, w in zip(Vv, Wd)), min(Vv) - mp.mpf(1) / 32)
    Rpt = mp.sqrt(8 * (S - mp.mpf(m - 1) / 32)) - mp.mpf(1) / 2
    print(f"\n== {name}: {m} cubics + {lin} lin, level {level}; S = {mp.nstr(S, 12)}, R(S) = {mp.nstr(Rpt, 8)}, max v = {mp.nstr(max(Vv), 6)}, "
          f"domain margin min(8v - (|w|+1/2)^2, v - 1/32) = {mp.nstr(dommarg, 4)}")
    print(f"   equalities (direct): E0 {mp.nstr(E0, 3)}, E1 {mp.nstr(E1, 3)}; odd parts {mp.nstr(O0, 3)}, {mp.nstr(O1, 3)}")
    # moment-form cross-check of E0, E1, PHI, and T_n^e via aud13 rows
    mom = {}
    for nm in ("E0", "E1", "PHI"):
        col, pc, c = AC.form(lin, nm)
        val = sum((mp.mpf(c0.numerator) / c0.denominator + mp.mpf(c1.numerator) / c1.denominator * S) *
                  sum(v ** a * (w * w) ** j for v, w in zip(Vv, Wd)) for (a, j), (c0, c1) in col.items())
        val += mp.mpf(pc.numerator) / pc.denominator * P + sum(mp.mpf(ci.numerator) / ci.denominator * S ** i for i, ci in enumerate(c))
        mom[nm] = val
    print(f"   moment forms vs direct: |E0m - E0| {mp.nstr(abs(mom['E0'] - E0), 3)}, |E1m - E1| {mp.nstr(abs(mom['E1'] - E1), 3)}, "
          f"|PHIm - Phi| {mp.nstr(abs(mom['PHI'] - Phi), 3)}")
    print(f"   CS: P - S M12/2 = {mp.nstr(P - S * M12 / 2, 6)}  ({'ok' if P <= S * M12 / 2 else 'VIOLATED'})")
    print(f"   PHI: Phi = {mp.nstr(Phi, 8)}  ({'ok' if Phi >= 0 else 'VIOLATED'});  (e2(beta)+e2(gamma))/2^10 = {mp.nstr(N3, 8)}, Phi - that = {mp.nstr(Phi - N3, 3)}")
    print(f"   N3L (all Bh, best Bh = Sb): Phi - 2^-10(Sb^2 - sum b^2 - sum om^2) = {mp.nstr(Phi - N3Lbest, 6)}  ({'ok' if Phi >= N3Lbest else 'VIOLATED'})")
    print(f"   N3U: 2^-10[Bbar Sb - sum b^2 + (Mtot-1) sum om^2] - Phi = {mp.nstr(N3Uval - Phi, 6)}  ({'ok' if N3Uval >= Phi else 'VIOLATED'})")
    print(f"   POS on xn t0 list: min Ve(t0) - 2^-8 t0^8 = {mp.nstr(pos[posmin_t0], 6)} at t0 = {posmin_t0}  ({'ok' if pos[posmin_t0] >= 0 else 'VIOLATED'}); "
          f"#violated t0 = {sum(1 for t in pos if pos[t] < 0)}; fine grid [-12, 12] min {mp.nstr(posgmin[0], 6)} at t = {mp.nstr(posgmin[1], 4)}")
    print(f"   Re(-7/5) = {mp.nstr(Re75, 8)};  S-identity: Re(-7/5) - lam Phi - RHS = {mp.nstr(Re75 - lamf * Phi - RHS, 3)};  "
          f"Phi < Phimax: {Phi < phimax};  RHS + lam Phimax = {mp.nstr(RHS + lamf * phimax, 6)}")
    # certificates containing S
    for c in certs[lin]:
        _, Ra, Rb, y, zl = c
        Sa, Sb_ = AC.SofR(Ra, m), AC.SofR(Rb, m)
        if not (Sa <= S <= Sb_):
            continue
        E = AC.eq_rows(lin)
        terms = [(yy, r) for yy, r in zip(y, E)]
        res_eq = [sum(1 for _ in Vv) - m, mp.mpf(0), E0, E1]       # mass, int v (S defined as sum v), E0, E1
        res_in = []
        for kind, prm, z in zl:
            terms.append((z, AC.ineq_row(lin, kind, prm)))
            if kind == "CS":
                res_in.append(P - S * M12 / 2)
            elif kind == "PHI":
                res_in.append(-Phi)
            elif kind == "N3L":
                Bh = mp.mpf(Fr(prm).numerator) / Fr(prm).denominator
                res_in.append(-(Phi - (2 * Bh * Sb - Bh ** 2 - sb2 - so2) / 1024))
            elif kind == "N3U":
                res_in.append(-(N3Uval - Phi))
            elif kind == "POS":
                t0 = Fr(prm)
                tt = mp.mpf(t0.numerator) / t0.denominator
                res_in.append(-(Ve(tt) - tt ** 8 / 256))
        col, pc, kap = AC.lin_comb(terms)

        def f2m(x):
            return mp.mpf(x.numerator) / x.denominator

        psi = sum((f2m(c0) + f2m(c1) * S) * sum(v ** a * (w * w) ** j for v, w in zip(Vv, Wd)) for (a, j), (c0, c1) in col.items())
        kS = f2m(kap[0]) + f2m(kap[1]) * S + f2m(kap[2]) * S ** 2
        lhs = psi + f2m(pc) * P - kS
        rhs = sum(f2m(yy) * r for yy, r in zip(y, res_eq)) + sum(f2m(z) * r for (_, _, z), r in zip(zl, res_in))
        viol = [(k, str(p), mp.nstr(r, 5)) for (k, p, z), r in zip(zl, res_in) if r > 0]
        print(f"   CERT R=[{Ra}, {Rb}] S=[{float(Sa):.5f}, {float(Sb_):.5f}] rows {[k + ('' if k in ('CS', 'PHI', 'N3U') else '@' + str(p)) for k, p, _ in zl]}: "
              f"Farkas identity lhs {mp.nstr(lhs, 10)} rhs {mp.nstr(rhs, 10)} (diff {mp.nstr(lhs - rhs, 3)}); sum_members Psi = {mp.nstr(psi, 6)} "
              f"(>= 0 expected), -kappa = {mp.nstr(-kS, 6)};  violated used rows: {viol if viol else 'NONE -> COVER REFUTED'}")


def main():
    certs = {0: [AC.parse_gp(l.strip()) for l in open("xn_cert9_l0.gp") if l.strip()],
             1: [AC.parse_gp(l.strip()) for l in open("xn_cert9_l1.gp") if l.strip()]}
    print("reconstructing m9 points")
    pts = sym_points()
    for p in pts:
        r = boundary_system(p[1], p[2], 0, 1)
        print(f"  {p[0]}: direct boundary residuals {[mp.nstr(x, 3) for x in r]}")
    pts.append(kcert_point("m9_kcert_9_0.out", "C90", 9, 0))
    pts.append(kcert_point("m9_kcert_9_1.out", "C81", 8, 1))
    pts.append(l2_point(9, 0))
    pts.append(l2_point(8, 1))
    for p in pts:
        evaluate(p, certs)


if __name__ == "__main__":
    main()
