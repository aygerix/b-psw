#!/usr/bin/env python3
"""aud13 (audit of candidate Theorem 13, k = 9): INDEPENDENT derivation of the reflection-even moment forms.

Method differs from the line (xn_k9mom.gp uses log/exp of Lam_n):
  T(s) = prod_i (1 + x_i(s)) mod s^9 with x_i = 1/(t g_i) at t = 1/s, x_i = O(s^3) for cubic members, so
  T_cubic = 1 + p1 + (p1^2 - p2)/2 (e_3(x) = O(s^9)), p_j = sum_i x_i^j, and T = T_cubic (1 + s/2)^lin.
  V(t) = sum_n T_n t^(8-n) is the polynomial part of H = (F - 1)/(t prod g).
Validation: V is recomputed by EXACT polynomial division (F - 1) div (t prod g) at random rational member sets,
both halves, both cases, and compared coefficient by coefficient with the moment forms.
Output: aud13_rows.json with, for lin = 0, 1, the even forms
  E0 = Ve(-1) - 2^-8, E1 = Ve'(-1) + 8 2^-8, PHI = Ve''(-1) - 56 2^-8, TE[n] = T_n^e (n = 0..8),
each as {col: {"a,j": [c0, c1]} (per-member v^a W^j, coefficient c0 + c1 S), pc: coefficient of P = M11^2/2, c: [c0, c1, c2]}.
Also checks the xn S-identity at t0 = -7/5 modulo span(E0, E1) and the S-bound arithmetic.
"""
import json
import random
import sys
from fractions import Fraction as Fr

import sympy as sp

K = 9
s, v, w, S = sp.symbols("s v w S")
R = sp.Rational
MS = {}


def M(a, b):
    key = (a, b)
    if key not in MS:
        MS[key] = sp.Symbol(f"M{a}_{b}")
    return MS[key]


def tomom(expr):
    """single-member polynomial in (v, w) -> linear form in moment symbols"""
    P = sp.Poly(sp.expand(expr), v, w)
    out = 0
    for (a, b), cf in P.terms():
        if a == 0:
            raise ValueError("v^0 term in a cubic-member quantity: " + str(expr))
        out += cf * M(a, b)
    return out


def member_series():
    p = R(3, 2) + w
    q = R(1, 2) + w + 2 * v
    d = [sp.Integer(1), -p]
    for n in range(2, 6):
        d.append(sp.expand(-p * d[n - 1] - q * d[n - 2]))
    x = {n + 3: sp.expand(v * d[n]) for n in range(6)}           # [s^n] x, n = 3..8
    x2 = {6: v ** 2, 7: sp.expand(2 * v ** 2 * d[1]), 8: sp.expand(v ** 2 * (2 * d[2] + d[1] ** 2))}
    return x, x2


def trunc(expr, n=K):
    P = sp.Poly(sp.expand(expr), s)
    return sum(cf * s ** e[0] for e, cf in P.terms() if e[0] < n)


def T_forms(lin):
    x, x2 = member_series()
    p1 = sum(tomom(x[n]) * s ** n for n in range(3, 9))
    p2 = sum(tomom(x2[n]) * s ** n for n in range(6, 9))
    Tc = trunc(1 + p1 + (p1 ** 2 - p2) / 2)
    T = trunc(Tc * (1 + s / 2) ** lin)
    P = sp.Poly(sp.expand(T), s)
    return [sp.expand(P.coeff_monomial(s ** n)) for n in range(K)]


def reflect(expr):
    sub = {sym: (-1) ** b * sym for (a, b), sym in MS.items()}
    return sp.expand(expr.xreplace(sub))


def even(expr):
    return sp.expand((expr + reflect(expr)) / 2)


def Vder_at(Tn, x0, d):
    """d-th derivative of V(t) = sum_n T_n t^(8-n) at t = x0"""
    tt = sp.Symbol("tt")
    V = sum(Tn[n] * tt ** (K - 1 - n) for n in range(K))
    return sp.expand(sp.diff(V, tt, d).subs(tt, x0))


def decompose(expr):
    """expr polynomial in moments -> (col, pc, c); M1_0 -> S"""
    e = sp.expand(expr.xreplace({M(1, 0): S}))
    syms = [sym for key, sym in MS.items() if key != (1, 0)]
    P = sp.Poly(e, *syms, S) if syms else sp.Poly(e, S)
    col, pc, c = {}, sp.Integer(0), [sp.Integer(0)] * 3
    gens = P.gens
    for mon, cf in P.terms():
        degs = dict(zip(gens, mon))
        sdeg = degs.get(S, 0)
        others = {g: d for g, d in degs.items() if g != S and d > 0}
        if not others:
            if sdeg > 2:
                raise ValueError("S degree > 2")
            c[sdeg] += cf
            continue
        if len(others) == 1:
            (g, d), = others.items()
            name = str(g)
            a, b = map(int, name[1:].split("_"))
            if (a, b) == (1, 1):
                if d != 2 or sdeg != 0:
                    raise ValueError("unexpected M1_1 term")
                pc += 2 * cf           # P = M11^2 / 2
                continue
            if d != 1 or sdeg > 1:
                raise ValueError(f"unexpected monomial {mon}")
            if b % 2:
                raise ValueError("odd moment in an even form")
            key = f"{a},{b // 2}"
            col.setdefault(key, [sp.Integer(0), sp.Integer(0)])
            col[key][sdeg] += cf
            continue
        raise ValueError(f"product of two non-S moments: {mon}")
    return col, pc, c


def ser(col, pc, c):
    return {"col": {k: [str(x) for x in vv] for k, vv in col.items() if any(x != 0 for x in vv)},
            "pc": str(pc), "c": [str(x) for x in c]}


# ---------------- exact direct check by polynomial division ----------------
def pmul(a, b):
    o = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            o[i + j] += x * y
    return o


def pdivmod(num, den):
    num = num[:]
    while num and num[-1] == 0:
        num.pop()
    q = [Fr(0)] * max(len(num) - len(den) + 1, 1)
    while len(num) >= len(den):
        f = num[-1] / den[-1]
        sh = len(num) - len(den)
        q[sh] = f
        for i, d in enumerate(den):
            num[sh + i] -= f * d
        num.pop()
        while num and num[-1] == 0:
            num.pop()
    return q, num


def direct_V(members, lin):
    """members: list of (a, c) rationals; returns V coefficients low -> high (exact)"""
    F = [Fr(1)]
    G = [Fr(0), Fr(1)]
    for a, c in members:
        F = pmul(F, [Fr(1), c + 2, a + c, a])
        G = pmul(G, [c + 2, a + c, a])
    for _ in range(lin):
        F = pmul(F, [Fr(1), Fr(2)])
        G = pmul(G, [Fr(2)])
    Fm1 = F[:]
    Fm1[0] -= 1
    q, r = pdivmod(Fm1, G)
    return q


def moments_of(members_vw):
    vals = {}
    for (a, b), sym in MS.items():
        vals[sym] = sum(Fr(vv) ** a * Fr(ww) ** b for vv, ww in members_vw)
    return vals


def validate(Tn_full, lin, ntrial=6, seed=1):
    rnd = random.Random(seed + lin)
    m = K - lin
    worst = 0
    for trial in range(ntrial):
        mem = []
        for _ in range(m):
            vv = Fr(rnd.randint(1, 400), rnd.randint(1, 97))
            uu = Fr(rnd.randint(-300, 300), rnd.randint(1, 89))
            mem.append((vv, uu))
        for refl in (0, 1):
            memr = [(vv, (1 - uu) if refl else uu) for vv, uu in mem]
            ac = [(1 / vv, uu / vv) for vv, uu in memr]
            Vq = direct_V(ac, lin)
            vals = moments_of([(vv, uu - Fr(1, 2)) for vv, uu in memr])
            for n in range(K):
                val = Tn_full[n].xreplace({sym: sp.Rational(x.numerator, x.denominator) for sym, x in vals.items()})
                val = Fr(int(sp.numer(val)), int(sp.denom(val)))
                d = Vq[K - 1 - n] if K - 1 - n < len(Vq) else Fr(0)
                if val != d:
                    print(f"  MISMATCH lin={lin} trial {trial} refl {refl} T_{n}: moment {float(val)} direct {float(d)}")
                    worst = 1
    print(f"  validation lin={lin}: {ntrial} random rational member sets x 2 halves, V coefficients "
          f"{'ALL EQUAL (exact)' if not worst else 'MISMATCH'}")
    return worst == 0


def main():
    out = {}
    allok = True
    for lin in (0, 1):
        Tn = T_forms(lin)
        allok &= validate(Tn, lin)
        Te = [even(x) for x in Tn]
        E0 = Vder_at(Te, -1, 0) - R(1, 256)
        E1 = Vder_at(Te, -1, 1) + R(8, 256)
        PHI = Vder_at(Te, -1, 2) - R(56, 256)
        forms = {"E0": E0, "E1": E1, "PHI": PHI}
        for n in range(K):
            forms[f"T{n}"] = Te[n]
        dec = {name: ser(*decompose(x)) for name, x in forms.items()}
        out[str(lin)] = dec
        print(f"lin = {lin}:")
        for name in ("E0", "E1", "PHI"):
            print(f"  {name} = {sp.expand(forms[name].xreplace({M(1, 0): S}))}")
        # S-identity check (xn): Re(-7/5) - lam Phi = RHS on span(E0, E1)
        t0 = R(-7, 5)
        Ve = sum(Te[n] * t0 ** (K - 1 - n) for n in range(K))
        Re = (Ve - t0 ** 8 / 256) / (t0 + 1) ** 2
        if lin == 0:
            lam, RHS = R(59, 110), -R(2, 55) * S ** 2 - R(788, 125) * S - R(2, 55) * M(2, 0) + R(39391023, 1100000)
            phimax, m = R(729, 32), 9
        else:
            lam, RHS = R(27, 50), -R(1, 25) * S ** 2 - R(678, 125) * S - R(1, 25) * M(2, 0) + R(1294023, 50000)
            phimax, m = R(585, 32), 8
        al, be = sp.symbols("al be")
        diff = sp.expand((Re - lam * PHI - al * E0 - be * E1).xreplace({M(1, 0): S}) - RHS)
        syms = [sym for key, sym in MS.items() if key != (1, 0)] + [S]
        eqs = sp.Poly(diff, *syms).coeffs()
        sol = sp.solve(eqs, [al, be], dict=True)
        print(f"  S-identity at t0 = -7/5 (lambda = {lam}): multipliers on (E0, E1) = {sol}  "
              f"{'HOLDS' if sol else 'FAILS'}")
        # e2(beta) bound check and S-bound: (c2) S^2 + (c1) S < c0 + lam*phimax with B0 >= 0
        P = sp.Poly(sp.expand(RHS.xreplace({M(2, 0): 0}) + lam * phimax), S)
        roots = [r for r in sp.solve(P.as_expr(), S) if r.is_real and r > 0]
        print(f"  S-bound: RHS + lam*Phimax > 0 requires S < {[sp.N(r, 12) for r in roots]}  (Phimax = {phimax})")
        Bbar = 18 * m + 2 * lin
        e2max = sp.binomial(m, 2) * 18 ** 2 + lin * m * 18 * 2
        print(f"  Phimax check: 2 e2max / 2^10 = {sp.Rational(2 * e2max, 1024)} vs {phimax}")
    with open("aud13_rows.json", "w") as f:
        json.dump(out, f, indent=1)
    print("validation:", "PASS" if allok else "FAIL")


if __name__ == "__main__":
    main()
