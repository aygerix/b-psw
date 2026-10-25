#!/usr/bin/env python3
"""xn line: RIGOROUS Farkas certificates for the k = 9 measure relaxation (reflection-even system + cuts).

Relaxation (see xn_lp9.py; all constraints are NECESSARY for a real cubic Williams family with k = 9):
  variables: a positive measure nu on D_S (counting measure of the m cubic members), and P = A1^2/2 >= 0;
  equalities: nu(D) = m, int v dnu = S, Es0 = 0, Es1 = 0;
  inequalities (flags): CS:  P <= S A2 / 2;  PHI: Phi >= 0;  N3L: tangent lower cuts;  N3U: upper cut;  POS: Ve(t0) >= 2^-8 t0^8.
Per-member column of every constraint is an affine function of S times the basis
  BASIS = [1, v, v^2, v^3, vW, v^2W, vW^2, 1/v, 1/v^2, W/v^2]   (W = w^2)
and every right-hand side is a polynomial in S of degree <= 2.
Certificate for an S-interval [Sa, Sb]: rational multipliers y (equalities), z >= 0 (inequalities) with
  (i)   Psi_S(W, v) = sum y_i col_i + sum z_j col_j >= 0 on D_{Sb} for S = Sa and S = Sb (Psi is affine in S, so this
        covers all S in [Sa, Sb]; D_S is contained in D_{Sb});
  (ii)  the P-column  sum y_i Pc_i + sum z_j Pc_j >= 0;
  (iii) kappa(S) = sum y_i rhs_i(S) + sum z_j rhs_j(S) < 0 for all S in [Sa, Sb].
Then (Farkas) no measure, hence no family, has S in [Sa, Sb].  (i) is verified EXACTLY: v = (r + 1/2)^2/8,
W = theta r^2, multiply by v^2 > 0, Bernstein expansion on boxes in (r, theta) with exact rationals, subdivision.
Usage: xn_farkas9.py lin Sa Sb flags [nv nw]     flags: subset of CS,PHI,N3L,N3U,POS joined by '+'
"""
import sys
from fractions import Fraction as Fr
from math import comb
import numpy as np
from scipy.optimize import linprog

BASIS = ["one", "v", "v2", "v3", "vW", "v2W", "vW2", "iv", "iv2", "Wiv2"]
NB = len(BASIS)


def model(lin):
    """returns m, K0, K1, KP, TE, T12 (exact Fractions)"""
    F = Fr
    if lin == 0:
        K0 = dict(SS=F(39, 8), SA2=F(1), SB0=F(-2), S=F(-321, 32), P=F(1), A2=F(-51, 4), B0=F(261, 8), A4=F(-9, 2), B2=F(57, 2), C0=F(-20), c=F(255, 256))
        K1 = dict(SS=F(-5, 2), S=F(351, 16), A2=F(12), B0=F(-27), A4=F(1), B2=F(-6), C0=F(4), c=F(-255, 32))
        KP = dict(SS=F(1), S=F(-209, 4), A2=F(-11), B0=F(23), c=F(1785, 32))
        TE = [dict(S=F(1)), dict(S=F(-3, 2)), dict(S=F(7, 4), A2=F(1), B0=F(-2)), dict(SS=F(1, 2), S=F(-15, 8), A2=F(-5, 2), B0=F(11, 2)),
              dict(SS=F(-3, 2), S=F(31, 16), A2=F(4), B0=F(-10), A4=F(1), B2=F(-6), C0=F(4)),
              dict(SS=F(23, 8), SA2=F(1), SB0=F(-2), S=F(-63, 32), P=F(1), A2=F(-21, 4), B0=F(121, 8), A4=F(-7, 2), B2=F(45, 2), C0=F(-16))]
        T12 = (F(0), F(0))
        m = 9
    else:
        K0 = dict(SS=F(31, 8), SA2=F(1), SB0=F(-2), S=F(-6), P=F(1), A2=F(-9), B0=F(191, 8), A4=F(-4), B2=F(51, 2), C0=F(-18), c=F(127, 256))
        K1 = dict(SS=F(-9, 4), S=F(15), A2=F(39, 4), B0=F(-89, 4), A4=F(1), B2=F(-6), C0=F(4), c=F(-143, 32))
        KP = dict(SS=F(1), S=F(-40), A2=F(-10), B0=F(21), c=F(1113, 32))
        TE = [dict(S=F(1)), dict(S=F(-1)), dict(S=F(1), A2=F(1), B0=F(-2)), dict(SS=F(1, 2), S=F(-1), A2=F(-2), B0=F(9, 2)),
              dict(SS=F(-5, 4), S=F(1), A2=F(11, 4), B0=F(-29, 4), A4=F(1), B2=F(-6), C0=F(4)),
              dict(SS=F(17, 8), SA2=F(1), SB0=F(-2), S=F(-1), P=F(1), A2=F(-13, 4), B0=F(81, 8), A4=F(-3), B2=F(39, 2), C0=F(-14))]
        T12 = (F(1, 2), F(0))
        m = 8
    return m, K0, K1, KP, TE, T12


MOM2BASIS = {"A2": "vW", "B0": "v2", "A4": "vW2", "B2": "v2W", "C0": "v3"}   # per-member value of each moment


class Row:
    """column = sum_b (c0[b] + c1[b] S) * basis_b ; P coefficient pc ; rhs = r0 + r1 S + r2 S^2"""

    def __init__(self):
        self.c0 = [Fr(0)] * NB
        self.c1 = [Fr(0)] * NB
        self.pc = Fr(0)
        self.rhs = [Fr(0)] * 3

    def add(self, other, f):
        for i in range(NB):
            self.c0[i] += f * other.c0[i]
            self.c1[i] += f * other.c1[i]
        self.pc += f * other.pc
        for i in range(3):
            self.rhs[i] += f * other.rhs[i]
        return self


def kform(K):
    """linear form of a moment expression K: returns Row with column part and 'constant' kept in rhs as NEGATED"""
    r = Row()
    for key, b in MOM2BASIS.items():
        if key in K:
            r.c0[BASIS.index(b)] += K[key]
    if "SA2" in K:
        r.c1[BASIS.index("vW")] += K["SA2"]
    if "SB0" in K:
        r.c1[BASIS.index("v2")] += K["SB0"]
    r.pc = K.get("P", Fr(0))
    # constant part c(S) = SS S^2 + S S + c  (stored as +c(S); callers move it)
    r.rhs = [K.get("c", Fr(0)), K.get("S", Fr(0)), K.get("SS", Fr(0))]
    return r


PARAMS = []
LAST_CERT = None


def build(lin, flags, Bh_list=None, t0_list=None):
    PARAMS.clear()
    m, K0, K1, KP, TE, T12 = model(lin)
    eqs, ineqs, names = [], [], []
    r = Row(); r.c0[0] = Fr(1); r.rhs = [Fr(m), Fr(0), Fr(0)]; eqs.append(r)
    r = Row(); r.c0[1] = Fr(1); r.rhs = [Fr(0), Fr(1), Fr(0)]; eqs.append(r)
    for K in (K0, K1):     # form(nu) + P pc + c(S) = 0  ->  form = -c(S)
        f = kform(K)
        f.rhs = [-x for x in f.rhs]
        eqs.append(f)
    phi = kform(KP)        # Phi = form + c(S)
    if "CS" in flags:      # P - (S/2) vW <= 0
        r = Row(); r.pc = Fr(1); r.c1[BASIS.index("vW")] = Fr(-1, 2); ineqs.append(r); names.append("CS"); PARAMS.append(("CS", 0))
    if "PHI" in flags:     # -Phi <= 0 : -form <= c(S)
        r = Row().add(phi, Fr(-1)); r.rhs = list(phi.rhs); ineqs.append(r); names.append("PHI"); PARAMS.append(("PHI", 0))
    Mtot = m + lin
    blin = 2 * lin
    # b = 2 + iv/2 ; b^2 = 4 + 2 iv + iv2/4 ; om^2 = Wiv2
    if "N3L" in flags:
        for Bh in Bh_list:
            # -Phi + 2^-10 [2 Bh sum b - Bh^2 - sum b^2 - sum om^2] <= 0
            r = Row().add(phi, Fr(-1))
            q = Fr(1, 1024)
            r.c0[0] += q * (2 * Bh * 2 - 4)
            r.c0[BASIS.index("iv")] += q * (2 * Bh * Fr(1, 2) - 2)
            r.c0[BASIS.index("iv2")] += q * (-Fr(1, 4))
            r.c0[BASIS.index("Wiv2")] += q * (-1)
            r.rhs = [phi.rhs[0] - q * (2 * Bh * blin - Bh * Bh - 4 * lin), phi.rhs[1], phi.rhs[2]]
            ineqs.append(r); names.append(f"N3L Bh={float(Bh):.3f}"); PARAMS.append(("N3L", Bh))
    if "N3U" in flags:
        Bbar = 18 * m + 2 * lin
        q = Fr(1, 512)
        r = Row().add(phi, Fr(1))
        # - q [ (Bbar b - b^2)/2 + (Mtot-1)/2 om^2 ]
        r.c0[0] -= q * (Bbar * 2 - 4) / 2
        r.c0[BASIS.index("iv")] -= q * (Bbar * Fr(1, 2) - 2) / 2
        r.c0[BASIS.index("iv2")] -= q * (-Fr(1, 4)) / 2
        r.c0[BASIS.index("Wiv2")] -= q * Fr(Mtot - 1, 2)
        r.rhs = [-phi.rhs[0] + q * (Bbar * blin - 4 * lin) / 2, -phi.rhs[1], -phi.rhs[2]]
        ineqs.append(r); names.append("N3U"); PARAMS.append(("N3U", 0))
    if "POS" in flags:
        TF = [kform(K) for K in TE]
        for t0 in t0_list:
            # -(sum_n Tn t0^(8-n)) <= (1 - 2^-8) t0^8 + T1 t0^7 + T2 t0^6 + sum_n cn(S) t0^(8-n)
            r = Row()
            rhs = [(1 - Fr(1, 256)) * t0 ** 8 + T12[0] * t0 ** 7 + T12[1] * t0 ** 6, Fr(0), Fr(0)]
            for j, f in enumerate(TF):
                p = t0 ** (8 - (j + 3))
                r.add(f, -p)
                for i in range(3):
                    rhs[i] += f.rhs[i] * p      # note: r.add also added -p * f.rhs; overwrite below
            r.rhs = rhs
            ineqs.append(r); names.append(f"POS t0={float(t0):.3f}"); PARAMS.append(("POS", t0))
    return m, eqs, ineqs, names


def basis_vals(W, v):
    return [1.0, v, v * v, v ** 3, v * W, v * v * W, v * W * W, 1 / v, 1 / v ** 2, W / v ** 2]


def col_matrix(rows, pts, S):
    """matrix (len(rows) x len(pts)) of column values at S (floats)"""
    B = np.array([basis_vals(W, v) for (W, v) in pts]).T     # NB x npts
    C0 = np.array([[float(x) for x in r.c0] for r in rows])
    C1 = np.array([[float(x) for x in r.c1] for r in rows])
    return (C0 + S * C1) @ B


def grid(m, Sb, nv, nw):
    vmax = Sb - (m - 1) / 32
    pts = []
    for v in np.linspace(1 / 32, vmax, nv):
        r = max(np.sqrt(8 * v) - 0.5, 0.0)
        for th in np.linspace(0, 1, nw):
            pts.append((th * r * r, v))
    return pts


def farkas_lp(lin, Sa, Sb, flags, nv=120, nw=9, Bh_list=None, t0_list=None, nSk=7, ybound=1e4, extra=None):
    m, eqs, ineqs, names = build(lin, flags, Bh_list, t0_list)
    pts = grid(m, Sb, nv, nw) + (extra or [])
    ne, ni = len(eqs), len(ineqs)
    nvar = ne + ni + 1          # y, z, eps
    A_ub, b_ub = [], []
    for S in (Sa, Sb):
        Me = col_matrix(eqs, pts, S)
        Mi = col_matrix(ineqs, pts, S) if ni else np.zeros((0, len(pts)))
        scale = 1 + np.abs(np.vstack([Me, Mi])).max(axis=0)
        for g in range(len(pts)):     # -(y.Me + z.Mi) + eps*scale <= 0
            row = np.zeros(nvar)
            row[:ne] = -Me[:, g]
            row[ne:ne + ni] = -Mi[:, g]
            row[-1] = scale[g]
            A_ub.append(row); b_ub.append(0.0)
    # P column >= 0
    row = np.zeros(nvar)
    row[:ne] = [-float(r.pc) for r in eqs]
    row[ne:ne + ni] = [-float(r.pc) for r in ineqs]
    A_ub.append(row); b_ub.append(-1e-6)
    # kappa(S) <= -1 on a grid of S
    for S in np.linspace(Sa, Sb, nSk):
        row = np.zeros(nvar)
        row[:ne] = [float(r.rhs[0] + r.rhs[1] * S + r.rhs[2] * S * S) for r in eqs]
        row[ne:ne + ni] = [float(r.rhs[0] + r.rhs[1] * S + r.rhs[2] * S * S) for r in ineqs]
        A_ub.append(row); b_ub.append(-1.0)
    c = np.zeros(nvar); c[-1] = -1
    bounds = [(-ybound, ybound)] * ne + [(0, ybound)] * ni + [(None, 1.0)]
    res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub), bounds=bounds, method="highs")
    return m, eqs, ineqs, names, res


# ---------------- exact verification ----------------
def poly_mul(a, b):
    out = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def poly_pow(a, e):
    out = [Fr(1)]
    for _ in range(e):
        out = poly_mul(out, a)
    return out


def psi_poly(coef):
    """coef: exact basis coefficients of Psi at fixed S.  Returns v^2 * Psi as dict theta_power -> poly in r."""
    vr = [Fr(1, 32), Fr(1, 8), Fr(1, 8)]          # v = (r + 1/2)^2 / 8 = 1/32 + r/8 + r^2/8
    Wr = [Fr(0), Fr(0), Fr(1)]                    # W = theta * r^2
    vp = [poly_pow(vr, e) for e in range(6)]
    # v^2 * basis:  one->v^2, v->v^3, v2->v^4, v3->v^5, vW->v^3 W, v2W->v^4 W, vW2->v^3 W^2, iv->v, iv2->1, Wiv2->W
    terms = {"one": (2, 0), "v": (3, 0), "v2": (4, 0), "v3": (5, 0), "vW": (3, 1), "v2W": (4, 1), "vW2": (3, 2),
             "iv": (1, 0), "iv2": (0, 0), "Wiv2": (0, 1)}
    out = {0: [Fr(0)], 1: [Fr(0)], 2: [Fr(0)]}
    for b, cf in zip(BASIS, coef):
        if cf == 0:
            continue
        ev, ew = terms[b]
        p = vp[ev]
        p = poly_mul(p, poly_pow(Wr, ew))
        p = [cf * x for x in p]
        cur = out[ew]
        if len(cur) < len(p):
            cur += [Fr(0)] * (len(p) - len(cur))
        for i, x in enumerate(p):
            cur[i] += x
    return out


def shift_scale(p, r0, h):
    """coefficients of p(r0 + h x) in x"""
    n = len(p)
    out = [Fr(0)] * n
    for i, a in enumerate(p):
        if a == 0:
            continue
        # a (r0 + h x)^i
        for j in range(i + 1):
            out[j] += a * comb(i, j) * r0 ** (i - j) * h ** j
    return out


def bernstein(p):
    n = len(p) - 1
    return [sum(Fr(comb(k, i), comb(n, i)) * p[i] for i in range(k + 1)) for k in range(n + 1)]


def verify_box(PP, r0, r1, t0=Fr(0), t1=Fr(1), depth=0, maxdepth=60, stats=None):
    """PP: dict j -> poly in r (theta^j coefficient, j <= 2).  Verify sum_j PP[j](r) theta^j > 0 on [r0,r1]x[t0,t1]."""
    h, ht = r1 - r0, t1 - t0
    deg = max(len(PP[j]) for j in PP) - 1
    # coefficient of y^i (theta = t0 + ht y) as a polynomial in r
    Q = {}
    for i in range(3):
        acc = [Fr(0)] * (deg + 1)
        for j in range(i, 3):
            f = comb(j, i) * t0 ** (j - i) * ht ** i
            for k, a in enumerate(PP[j]):
                acc[k] += f * a
        Q[i] = bernstein(shift_scale(acc, r0, h))
    mn = None
    for k in range(deg + 1):
        for l in range(3):
            val = sum(Fr(comb(l, i), comb(2, i)) * Q[i][k] for i in range(l + 1))
            mn = val if mn is None or val < mn else mn
    if stats is not None:
        stats["boxes"] += 1
        if stats["boxes"] > 200000:
            return False
    if mn > 0:
        return True
    if depth >= maxdepth:
        return False
    if depth % 3 == 2:      # split theta every third level, r otherwise
        tm = (t0 + t1) / 2
        return (verify_box(PP, r0, r1, t0, tm, depth + 1, maxdepth, stats) and
                verify_box(PP, r0, r1, tm, t1, depth + 1, maxdepth, stats))
    rm = (r0 + r1) / 2
    return (verify_box(PP, r0, rm, t0, t1, depth + 1, maxdepth, stats) and
            verify_box(PP, rm, r1, t0, t1, depth + 1, maxdepth, stats))


def rat(x, D=10 ** 6):
    return Fr(x).limit_denominator(D)


def S_of_R(R, m):
    return (R + Fr(1, 2)) ** 2 / 8 + Fr(m - 1, 32)


def R_of_S(S, m, up):
    x = np.sqrt(8 * (S - (m - 1) / 32)) - 0.5
    q = Fr(int(np.ceil(x * 4096)) if up else int(np.floor(x * 4096)), 4096)
    return max(q, Fr(0))


def certify(lin, Sa, Sb, flags, nv=120, nw=9, nBh=40, nt0=36, ybound=1e4):
    m = 9 - lin
    return certify_R(lin, R_of_S(Sa, m, False), R_of_S(Sb, m, True), flags, nv, nw, nBh, nt0, ybound)


def certify_R(lin, Ra, Rb, flags, nv=120, nw=9, nBh=40, nt0=36, ybound=1e4):
    m = 9 - lin
    Sa, Sb = S_of_R(Ra, m), S_of_R(Rb, m)          # exact rational endpoints; the cap v <= Sb - (m-1)/32 <-> r <= Rb
    lo = 2 * m + Fr(m * m) / (2 * Sb) + 2 * lin
    Bh_list = [lo + (18 * m + 2 * lin - lo) * Fr(i, nBh - 1) for i in range(nBh)] if "N3L" in flags else None
    t0_list = [Fr(-4) + Fr(7 * i, nt0 - 1) for i in range(nt0)] if "POS" in flags else None
    extra = []
    for rnd in range(6):
        ok, msg, bad = certify_once(lin, Sa, Sb, Ra, Rb, flags, nv, nw, Bh_list, t0_list, ybound, extra)
        if ok or bad is None:
            return ok, msg
        extra += bad
    return ok, msg + f" (after {rnd + 1} cutting-plane rounds)"


def certify_once(lin, Sa, Sb, Ra, Rb, flags, nv, nw, Bh_list, t0_list, ybound, extra):
    m = 9 - lin
    m, eqs, ineqs, names, res = farkas_lp(lin, float(Sa), float(Sb), flags, nv, nw, Bh_list, t0_list, ybound=ybound, extra=extra)
    if res.status != 0:
        return False, f"Farkas LP status {res.status}: {res.message}", None
    eps = res.x[-1]
    if eps <= 0:
        return False, f"no margin (eps = {eps:.3e})", None
    ne, ni = len(eqs), len(ineqs)
    y = [rat(x) for x in res.x[:ne]]
    z = [max(Fr(0), rat(x)) for x in res.x[ne:ne + ni]]
    # (ii) P column
    pcol = sum(yi * r.pc for yi, r in zip(y, eqs)) + sum(zi * r.pc for zi, r in zip(z, ineqs))
    if pcol < 0:
        return False, f"P column negative after rounding: {float(pcol)}", None
    # (iii) kappa(S) < 0 on [Sa, Sb]
    k0 = sum(yi * r.rhs[0] for yi, r in zip(y, eqs)) + sum(zi * r.rhs[0] for zi, r in zip(z, ineqs))
    k1 = sum(yi * r.rhs[1] for yi, r in zip(y, eqs)) + sum(zi * r.rhs[1] for zi, r in zip(z, ineqs))
    k2 = sum(yi * r.rhs[2] for yi, r in zip(y, eqs)) + sum(zi * r.rhs[2] for zi, r in zip(z, ineqs))
    cand = [Sa, Sb]
    if k2 != 0:
        sv = -k1 / (2 * k2)
        if Sa < sv < Sb:
            cand.append(sv)
    kmax = max(k0 + k1 * s + k2 * s * s for s in cand)
    if kmax >= 0:
        return False, f"kappa not negative: max {float(kmax)}", None
    # (i) Psi >= 0 on D_{Sb} at S = Sa and S = Sb
    R = Rb
    stats = {"boxes": 0}
    for S in (Sa, Sb):
        coef = [sum(yi * (r.c0[b] + r.c1[b] * S) for yi, r in zip(y, eqs)) + sum(zi * (r.c0[b] + r.c1[b] * S) for zi, r in zip(z, ineqs))
                for b in range(NB)]
        PP = psi_poly(coef)
        if not verify_box(PP, Fr(0), R, stats=stats):
            # locate the bad region numerically and return points for the next cutting-plane round
            Rf = float(R)
            pf = {j: [float(x) for x in PP[j]] for j in PP}
            bad = []
            rs = np.linspace(0, Rf, 2001)
            for th in np.linspace(0, 1, 41):
                vals = sum(np.polyval(pf[j][::-1], rs) * th ** j for j in pf)
                for i in np.argsort(vals)[:3]:
                    rr = rs[i]
                    v = (rr + 0.5) ** 2 / 8
                    bad.append((th * rr * rr, v))
            return False, f"Psi positivity not verified at S = {float(S)} (boxes {stats['boxes']})", bad
    used = [names[j] for j in range(ni) if z[j] > 0]
    global LAST_CERT
    LAST_CERT = dict(lin=lin, Ra=Ra, Rb=Rb, y=y, z=[(PARAMS[j][0], PARAMS[j][1], z[j]) for j in range(ni) if z[j] > 0])
    return True, f"CERTIFIED S in [{Sa}, {Sb}] = [{float(Sa):.6f}, {float(Sb):.6f}]: eps {eps:.3e}, kappa max {float(kmax):.4g}, boxes {stats['boxes']}, active ineqs {len(used)}: {used[:6]}", None


if __name__ == "__main__":
    lin = int(sys.argv[1])
    Sa, Sb = float(sys.argv[2]), float(sys.argv[3])
    flags = sys.argv[4].split("+") if len(sys.argv) > 4 and sys.argv[4] != "none" else []
    nv = int(sys.argv[5]) if len(sys.argv) > 5 else 120
    nw = int(sys.argv[6]) if len(sys.argv) > 6 else 9
    ok, msg = certify(lin, Sa, Sb, flags, nv, nw)
    print(f"lin={lin} S in [{Sa}, {Sb}] flags={flags}: {msg}", flush=True)
