#!/usr/bin/env python3
"""m9 line: numerical feasibility of the four boundary equations (real cubic Williams families, k = 9, 11).

Unknowns per cubic member: (rho, z), r = exp(rho) > 0, v = (r + 1/2)^2 / 8, u = 1/2 + r tanh(z).
This parametrises EXACTLY the admissible domain (|u - 1/2| + 1/2)^2 < 8 v (which implies v > 1/32).
Equations: T(-1) = 2^(1-k), T'(-1) = 0 for the family and for its reflection u -> 1 - u, where
T(s) = prod_i (1 + x_i(s)) * (1 + s/2)^lin mod s^k, x = v s^3 / (1 + (1+u) s + (u + 2v) s^2).
(Equivalent to V(-1) = 2^(1-k), V'(-1) = (1-k) 2^(1-k); identity checked exactly in m9_check.gp.)
Usage: m9_feas.py k lin nstarts seed
"""
import sys
import numpy as np
from scipy.optimize import least_squares


def xcoef(u, v, k):
    """coefficients 0..k-1 of x(s) mod s^k"""
    c = np.zeros(k)
    # series of 1/(1 + (1+u)s + (u+2v)s^2)
    d = np.zeros(k)
    d[0] = 1.0
    if k > 1:
        d[1] = -(1 + u)
    for j in range(2, k):
        d[j] = -(1 + u) * d[j - 1] - (u + 2 * v) * d[j - 2]
    for n in range(3, k):
        c[n] = v * d[n - 3]
    return c


def mulser(a, b, k):
    return np.convolve(a, b)[:k]


def Tvals(U, W, k, lin):
    T = np.zeros(k)
    T[0] = 1.0
    for u, v in zip(U, W):
        f = xcoef(u, v, k)
        f[0] += 1.0
        T = mulser(T, f, k)
    if lin:
        T = mulser(T, np.array([1.0, 0.5] + [0.0] * (k - 2)), k)
    n = np.arange(k)
    sg = (-1.0) ** n
    return np.dot(T, sg), np.dot(n[1:] * T[1:], sg[:-1])  # T(-1), T'(-1)


def unpack(p):
    m = len(p) // 2
    rho, z = p[:m], p[m:]
    r = np.exp(rho)
    v = (r + 0.5) ** 2 / 8
    u = 0.5 + r * np.tanh(z)
    return u, v


def resid(p, k, lin):
    u, v = unpack(p)
    c0 = 2.0 ** (1 - k)
    a0, a1 = Tvals(u, v, k, lin)
    b0, b1 = Tvals(1 - u, v, k, lin)
    return np.array([(a0 - c0) / c0, a1 / c0, (b0 - c0) / c0, b1 / c0])


def margins(u, v):
    """interior margins: sqrt(8v) - (|u-1/2| + 1/2) > 0 ; v - 1/32 > 0 ; min pairwise distance"""
    mg = np.sqrt(8 * v) - (np.abs(u - 0.5) + 0.5)
    pts = np.stack([u, v], 1)
    dmin = min(np.linalg.norm(pts[i] - pts[j]) for i in range(len(u)) for j in range(i + 1, len(u)))
    return mg.min(), (v - 1 / 32).min(), dmin


def main():
    k, lin, ns, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    m = k - lin
    rng = np.random.default_rng(seed)
    nok = 0
    for it in range(ns):
        p0 = np.concatenate([rng.normal(0, 1.5, m), rng.normal(0, 1.0, m)])
        try:
            sol = least_squares(resid, p0, args=(k, lin), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        except Exception as e:  # overflow etc.
            continue
        f = np.abs(sol.fun).max()
        u, v = unpack(sol.x)
        mg, mv, dm = margins(u, v)
        if f < 1e-10:
            nok += 1
            J = sol.jac
            sv = np.linalg.svd(J, compute_uv=False)
            print(f"start {it}: maxres {f:.2e} margin {mg:.4f} v-1/32 {mv:.4f} dmin {dm:.4f} "
                  f"Jac sv min {sv.min():.3e} rank {np.sum(sv > 1e-8 * sv.max())}")
            print("   u =", np.array2string(u, precision=12, max_line_width=250))
            print("   v =", np.array2string(v, precision=12, max_line_width=250))
            print("   S =", v.sum(), " a = 1/v:", np.array2string(1 / v, precision=6, max_line_width=250))
        else:
            print(f"start {it}: fail maxres {f:.2e} (margin {mg:.3e})")
    print("converged:", nok, "of", ns)


if __name__ == "__main__":
    main()
