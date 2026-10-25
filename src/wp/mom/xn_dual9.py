#!/usr/bin/env python3
"""xn line: Lagrangian (single-atom) dual bound for k = 9, 9 cubics, reflection-EVEN system.

Necessary for a family (moment form from xn_k9mom.gp, w = u - 1/2, nu' = sum v_i delta_(w_i, v_i), mass S):
  Es0 := (E0 + E0r)/2 = 0,   Es1 := (E1 + E1r)/2 = 0,   Phi := (D2 + D2r)/2 = 2^(-k-1)(e2(beta) + e2(gamma)) > 0.
Dual: for alpha <= 0, beta real, at fixed S:
  Phi = Phi + alpha Es0 + beta Es1 = int psi_S dnu' + c(S) + alpha A1^2/2  <=  S * sup_{D_S} psi_S + c(S),
  D_S = {|w| < sqrt(8v) - 1/2, 1/32 < v <= S - (m-1)/32}.
If  min_{alpha<=0, beta} [S sup psi_S + c(S)] < 0  for every S, the even system + Phi > 0 is infeasible.
Numerical sup (grid + refinement): COMPUTED only; a rigorous version needs a certified sup.
Usage: xn_dual9.py Smin Smax nS
"""
import sys
import numpy as np
from scipy.optimize import minimize

m = 9


def psi(w, v, S, al, be):
    w2 = w * w
    return (w2 * w2 * (be - 4.5 * al) + w2 * (-11 - 51 / 4 * al + 12 * be + al * S) + v * w2 * (57 / 2 * al - 6 * be)
            + v * (23 + 261 / 8 * al - 27 * be - 2 * al * S) + v * v * (4 * be - 20 * al))


def cS(S, al, be):
    return (S * S - 209 / 4 * S + 1785 / 32 + al * (39 / 8 * S * S - 321 / 32 * S + 255 / 256)
            + be * (-5 / 2 * S * S + 351 / 16 * S - 255 / 32))


def grid(S, n=400):
    vmax = S - (m - 1) / 32
    vs = np.linspace(1 / 32, vmax, n)
    pts = []
    for v in vs:
        r = max(np.sqrt(8 * v) - 0.5, 0)
        ws = np.linspace(0, r, 60)     # psi is even in w
        pts.append(np.stack([ws, np.full_like(ws, v)], 1))
    return np.concatenate(pts)


def bound(x, S, G):
    al, be = -np.exp(x[0]) if x[0] < 50 else -1e20, x[1]
    vals = psi(G[:, 0], G[:, 1], S, al, be)
    return S * vals.max() + cS(S, al, be)


def bound_al(al, be, S, G):
    vals = psi(G[:, 0], G[:, 1], S, al, be)
    return S * vals.max() + cS(S, al, be)


Smin, Smax, nS = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
for S in np.linspace(Smin, Smax, nS):
    if S <= m / 32:
        continue
    G = grid(S)
    best = (np.inf, None)
    for al0 in [0.0, -0.01, -0.1, -1.0, -5.0]:
        for be0 in [-5.0, -1.0, 0.0, 1.0, 5.0]:
            f = lambda x: bound_al(min(x[0], 0.0), x[1], S, G)
            r = minimize(f, [al0, be0], method="Nelder-Mead", options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 4000})
            if r.fun < best[0]:
                best = (r.fun, (min(r.x[0], 0.0), r.x[1]))
    print(f"S = {S:.4f}: dual bound {best[0]:.6g} at alpha = {best[1][0]:.6g}, beta = {best[1][1]:.6g}", flush=True)
