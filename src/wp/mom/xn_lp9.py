#!/usr/bin/env python3
"""xn line: LP (measure) relaxation for k = 9 of  {reflection-even boundary equations, Phi_even >= 0}.

Moment forms (xn_k9mom.gp; w = u - 1/2, sums over the m cubic members; nu = counting measure):
  S = sum v, A1 = sum v w, A2 = sum v w^2, A4 = sum v w^4, B0 = sum v^2, B2 = sum v^2 w^2, C0 = sum v^3, P = A1^2/2.
  lin = 0 (m = 9):
    Es0 = 39/8 S^2 + S (A2 - 2 B0 - 321/32) + P - 51/4 A2 + 261/8 B0 - 9/2 A4 + 57/2 B2 - 20 C0 + 255/256 = 0
    Es1 = -5/2 S^2 + 351/16 S + 12 A2 - 27 B0 + A4 - 6 B2 + 4 C0 - 255/32 = 0
    Phi = S^2 - 209/4 S - 11 A2 + 23 B0 + 1785/32     [= (D2 + D2r)/2 = 2^(-10)(e2(beta) + e2(gamma)) > 0]
  lin = 1 (m = 8, plus 2t + 1):
    Es0 = 31/8 S^2 + S (A2 - 2 B0 - 6) + P - 9 A2 + 191/8 B0 - 4 A4 + 51/2 B2 - 18 C0 + 127/256 = 0
    Es1 = -9/4 S^2 + 15 S + 39/4 A2 - 89/4 B0 + A4 - 6 B2 + 4 C0 - 143/32 = 0
    Phi = S^2 - 40 S - 10 A2 + 21 B0 + 1113/32
Relaxation at fixed S: nu any positive measure on D_S = {0 <= w < sqrt(8v) - 1/2, 1/32 < v <= S - (m-1)/32}
(even system: only w^2 enters) with nu(D_S) = m and int v dnu = S; P free in [0, S A2 / 2] (Cauchy-Schwarz).
LP maximises Phi.  max < 0 at this S  =>  no family with this S (up to grid discretisation: COMPUTED).
Optional N3 cuts (flag n3 = 1): the exact second-order equations D2 = 2^-k e2(beta), D2r = 2^-k e2(gamma) give
  Phi = 2^-10 (e2(beta) + e2(gamma)) = 2^-9 (e2(b) + e2(om)),  b_i = 2 + 1/(2 v_i), om_i = w_i / v_i
      >= 2^-10 [ 2 Bh sum b - Bh^2 - sum b^2 - sum om^2 ]   for EVERY real Bh  (tangent cut, (sum om)^2 >= 0 dropped).
  (lin = 1: the linear member has beta = gamma = 2, i.e. b = 2, om = 0, and k = 9 throughout.)
Optional positivity cuts (flag pos = 1): Ve(t0) - 2^-8 t0^8 >= 0 for t0 on a grid in [-4, 3], where
  Ve(t) = sum_n Tn^e t^(8-n) is the reflection-average of Vpoly (each of V, V~ is >= 2^(1-k) t^(k-1) on R, PROVED:
  V = 2^(1-k) sum_j t^(k-1-2j) (t+1)^(2j) e_2j(h/g) with h_i, g_i > 0).  Tn^e from xn_k9T.gp.
Optional N3 UPPER cut (flag n3u = 1): with a_i = 1/v_i < 32, b_i = 2 + a_i/2 <= 18, om_i = w_i a_i,
  Phi = 2^-9 (e2(b) + e2(om)) <= 2^-9 [ (Bbar sum b - sum b^2)/2 + (M - 1)/2 sum om^2 ],  Bbar = 18 M + 2 lin,
  M = number of members carrying b, om (cubics + linear); (sum b)^2 <= Bbar sum b, e2(om) <= (M-1)/2 sum om^2.
Usage: xn_lp9.py lin Smin Smax nS [nv nw n3 pos n3u]
"""
import sys
import numpy as np
from scipy.optimize import linprog

lin = int(sys.argv[1])
Smin, Smax, nS = float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
nv = int(sys.argv[5]) if len(sys.argv) > 5 else 300
nw = int(sys.argv[6]) if len(sys.argv) > 6 else 40
n3 = int(sys.argv[7]) if len(sys.argv) > 7 else 0
pos = int(sys.argv[8]) if len(sys.argv) > 8 else 0
n3u = int(sys.argv[9]) if len(sys.argv) > 9 else 0
sense = sys.argv[10] if len(sys.argv) > 10 else "max"     # "min": minimise Phi instead
if lin == 0:   # even parts T3e..T8e: dict of (SS, S, A2, B0, A4, B2, C0, SA2, SB0, P, c)
    TE = [dict(S=1, c=0), dict(S=-3 / 2, c=0), dict(S=7 / 4, A2=1, B0=-2, c=0), dict(SS=1 / 2, S=-15 / 8, A2=-5 / 2, B0=11 / 2, c=0),
          dict(SS=-3 / 2, S=31 / 16, A2=4, B0=-10, A4=1, B2=-6, C0=4, c=0),
          dict(SS=23 / 8, SA2=1, SB0=-2, S=-63 / 32, P=1, A2=-21 / 4, B0=121 / 8, A4=-7 / 2, B2=45 / 2, C0=-16, c=0)]
    T12 = (0.0, 0.0)
else:
    TE = [dict(S=1, c=0), dict(S=-1, c=0), dict(S=1, A2=1, B0=-2, c=0), dict(SS=1 / 2, S=-1, A2=-2, B0=9 / 2, c=0),
          dict(SS=-5 / 4, S=1, A2=11 / 4, B0=-29 / 4, A4=1, B2=-6, C0=4, c=0),
          dict(SS=17 / 8, SA2=1, SB0=-2, S=-1, P=1, A2=-13 / 4, B0=81 / 8, A4=-3, B2=39 / 2, C0=-14, c=0)]
    T12 = (0.5, 0.0)
m = 9 - lin
if lin == 0:
    K0 = dict(SS=39 / 8, SA2=1, SB0=-2, S=-321 / 32, P=1, A2=-51 / 4, B0=261 / 8, A4=-9 / 2, B2=57 / 2, C0=-20, c=255 / 256)
    K1 = dict(SS=-5 / 2, S=351 / 16, A2=12, B0=-27, A4=1, B2=-6, C0=4, c=-255 / 32)
    KP = dict(SS=1, S=-209 / 4, A2=-11, B0=23, c=1785 / 32)
else:
    K0 = dict(SS=31 / 8, SA2=1, SB0=-2, S=-6, P=1, A2=-9, B0=191 / 8, A4=-4, B2=51 / 2, C0=-18, c=127 / 256)
    K1 = dict(SS=-9 / 4, S=15, A2=39 / 4, B0=-89 / 4, A4=1, B2=-6, C0=4, c=-143 / 32)
    KP = dict(SS=1, S=-40, A2=-10, B0=21, c=1113 / 32)


def pts(S):
    vmax = S - (m - 1) / 32
    out = []
    for v in np.linspace(1 / 32 + 1e-9, vmax, nv):
        r = max(np.sqrt(8 * v) - 0.5, 0.0)
        for w in np.linspace(0, r, nw):
            out.append((w, v))
    return np.array(out)


def rowfuncs(G):
    w, v = G[:, 0], G[:, 1]
    return dict(one=np.ones_like(v), S=v, A2=v * w * w, A4=v * w ** 4, B0=v * v, B2=v * v * w * w, C0=v ** 3)


def lin_row(K, F, S):
    """coefficients on nu (vector over grid) of the S-dependent linear part, and the constant"""
    row = np.zeros(len(F["S"]))
    for key in ("A2", "B0", "A4", "B2", "C0"):
        if key in K:
            row += K[key] * F[key]
    if "SA2" in K:
        row += K["SA2"] * S * F["A2"]
    if "SB0" in K:
        row += K["SB0"] * S * F["B0"]
    const = K.get("SS", 0) * S * S + K.get("S", 0) * S + K["c"]
    return row, const


for S in np.linspace(Smin, Smax, nS):
    if S <= m / 32 + 1e-6:
        continue
    G = pts(S)
    F = rowfuncs(G)
    r0, c0 = lin_row(K0, F, S)
    r1, c1 = lin_row(K1, F, S)
    rp, cp = lin_row(KP, F, S)
    n = len(G)
    # variables: nu (n), P (1)
    A_eq = np.zeros((4, n + 1))
    b_eq = np.zeros(4)
    A_eq[0, :n] = F["one"]; b_eq[0] = m
    A_eq[1, :n] = F["S"]; b_eq[1] = S
    A_eq[2, :n] = r0; A_eq[2, n] = K0["P"]; b_eq[2] = -c0
    A_eq[3, :n] = r1; b_eq[3] = -c1
    # P - S A2 / 2 <= 0
    rows, rhs = [], []
    r = np.zeros(n + 1); r[:n] = -S * F["A2"] / 2; r[n] = 1
    rows.append(r); rhs.append(0.0)
    if n3:
        w_, v_ = G[:, 0], G[:, 1]
        bb = 2 + 1 / (2 * v_)
        om2 = (w_ / v_) ** 2
        blin = 2.0 * lin       # linear member: b = 2 (constant contribution to sum b), b^2 = 4
        lo = 2 * m + m * m / (2 * S) + blin
        for Bh in np.linspace(lo, 18 * m + blin, 80):
            # cut: -Phi + 2^-10 [2 Bh (sum b) - Bh^2 - sum b^2 - sum om^2] <= 0
            r = np.zeros(n + 1)
            r[:n] = -rp + 2.0 ** -10 * (2 * Bh * bb - bb * bb - om2)
            rows.append(r)
            rhs.append(cp - 2.0 ** -10 * (2 * Bh * blin - Bh * Bh - 4.0 * lin))
    if n3u:
        w_, v_ = G[:, 0], G[:, 1]
        bb = 2 + 1 / (2 * v_)
        om2 = (w_ / v_) ** 2
        Mtot = m + lin
        Bbar = 18 * m + 2 * lin
        # Phi - 2^-9 [ (Bbar (sum b) - sum b^2)/2 + (Mtot-1)/2 sum om^2 ] <= 0   (linear member: b = 2, om = 0)
        r = np.zeros(n + 1)
        r[:n] = rp - 2.0 ** -9 * ((Bbar * bb - bb * bb) / 2 + (Mtot - 1) / 2 * om2)
        rows.append(r)
        rhs.append(-cp + 2.0 ** -9 * (Bbar * 2 * lin - 4 * lin) / 2)
    if pos:
        TR = [lin_row(K, F, S) for K in TE]
        for t0 in np.linspace(-4, 3, 71):
            # -(sum_n Tn t0^(8-n) + (1 - 2^-8) t0^8) <= 0 ; n = 3..8 from TE, n = 1, 2 constants T12
            r = np.zeros(n + 1)
            cst = (1 - 2.0 ** -8) * t0 ** 8 + T12[0] * t0 ** 7 + T12[1] * t0 ** 6
            for j, (rw, cc) in enumerate(TR):
                nn = j + 3
                r[:n] -= rw * t0 ** (8 - nn)
                cst += cc * t0 ** (8 - nn)
                if "P" in TE[j]:
                    r[n] -= TE[j]["P"] * t0 ** (8 - nn)
            rows.append(r)
            rhs.append(cst)
    A_ub = np.array(rows)
    b_ub = np.array(rhs)
    cvec = np.zeros(n + 1); cvec[:n] = -rp if sense == "max" else rp
    res = linprog(cvec, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=[(0, None)] * (n + 1), method="highs")
    if res.status == 2:
        print(f"S = {S:.4f}: LP INFEASIBLE (even system has no measure solution)", flush=True)
    elif res.status == 0:
        val = (-res.fun if sense == "max" else res.fun) + cp
        sup = np.nonzero(res.x[:n] > 1e-9)[0]
        desc = "; ".join(f"({G[i,0]:.3f},{G[i,1]:.4f})x{res.x[i]:.3f}" for i in sup[:6])
        if len(sys.argv) > 11 and sys.argv[11] == "dual":
            el, il = res.eqlin.marginals, res.ineqlin.marginals
            print("   eq duals [count, S, Es0, Es1]:", np.array2string(el, precision=6))
            labels = ["CS(P)"] + ([f"N3cut{j}" for j in range(80)] if n3 else []) + (["N3upper"] if n3u else []) + ([f"pos t0={t0:.2f}" for t0 in np.linspace(-4, 3, 71)] if pos else [])
            for j in np.nonzero(np.abs(il) > 1e-9)[0]:
                print(f"   ineq dual {labels[j]}: {il[j]:.6g}")
        print(f"S = {S:.4f}: LP {sense} Phi = {val:.6g}   support[{len(sup)}]: {desc}  P={res.x[n]:.4g}", flush=True)
    else:
        print(f"S = {S:.4f}: LP status {res.status} {res.message}", flush=True)
