#!/usr/bin/env python3
"""m9 line: the Taylor hierarchy beyond the four boundary equations (real cubic Williams families).

For members f_j = (v_j + t Q_j(t)) / v_j  (cubic: Q = t^2 + (1+u) t + u + 2v; linear 2t+1: Q = 1, v = 1/2),
H := (F - 1) / (t prod_j (Q_j / v_j)) = t^(k-1) [prod(1 + x_j) - prod x_j] is a rational function with polynomial part V.
The "-" half of the Williams condition  <=>  H = V  <=>  R := (F-1)/t mod prod g_j = 0 (deg R < 2m)
  <=>  the Taylor coefficients r_0, ..., r_(2m-1) of H - V at t = -1 vanish.
r_0 = 0 and r_1 = 0 are exactly the two boundary equations (V(-1), V'(-1)).  The "+" half is the same for the
reflected family u -> 1 - u.  Level J of the hierarchy: r_0..r_J = 0 for both halves (2(J+1) equations).
This script tests REAL feasibility of level J in the exact admissible domain, by least squares from random starts.
Usage: m9_hier.py k lin J nstarts seed [check]
"""
import sys
import numpy as np
from numpy.polynomial import polynomial as P
from scipy.optimize import least_squares
sys.path.insert(0, '.')
from m9_feas import Tvals


def Vcoef(U, W, k, lin):
    """coefficients of V(t) (ascending powers) from the Laurent part"""
    T = np.zeros(k)
    T[0] = 1.0
    for u, v in zip(U, W):
        d = np.zeros(k)
        d[0] = 1.0
        if k > 1:
            d[1] = -(1 + u)
        for j in range(2, k):
            d[j] = -(1 + u) * d[j - 1] - (u + 2 * v) * d[j - 2]
        f = np.zeros(k)
        f[0] = 1.0
        f[3:] = v * d[:k - 3]
        T = np.convolve(T, f)[:k]
    if lin:
        T = np.convolve(T, np.array([1.0, 0.5] + [0.0] * (k - 2)))[:k]
    # V(t) = sum_j T_j t^(k-1-j)  -> ascending coefficient of t^i is T_(k-1-i)
    return T[::-1].copy()


def shift(c, a):
    """coefficients of p(a + z) in z, p ascending coefficients c"""
    n = len(c)
    out = np.zeros(n)
    # Horner in z: p(a+z)
    for co in c[::-1]:
        out = np.convolve(out, np.array([a, 1.0]))[:n]
        out[0] += co
    return out


def taylor_resid(U, W, k, lin, N):
    """r_0..r_(N-1): Taylor coefficients of H - V at t = -1, z = t + 1.
    x_j(-1+z) = v / D(z), D = z^3 + (u-2) z^2 + (2v-u+1) z - 2v  (cubic);  x_0 = 1/(2(z-1))  (linear member).
    H = (z-1)^(k-1) [prod(1+x) - prod x];  V(-1+z) by Taylor shift of the Laurent polynomial part."""
    def inv(dc):
        out = np.zeros(N)
        out[0] = 1.0 / dc[0]
        for n in range(1, N):
            acc = 0.0
            for l in range(1, min(n, len(dc) - 1) + 1):
                acc += dc[l] * out[n - l]
            out[n] = -acc / dc[0]
        return out
    xs = []
    for u, v in zip(U, W):
        xs.append(v * inv(np.array([-2 * v, 2 * v - u + 1, u - 2, 1.0])))
    if lin:
        xs.append(0.5 * inv(np.array([-1.0, 1.0])))
    A = np.zeros(N); A[0] = 1.0
    B = np.zeros(N); B[0] = 1.0
    for x in xs:
        one = x.copy(); one[0] += 1.0
        A = np.convolve(A, one)[:N]
        B = np.convolve(B, x)[:N]
    G = A - B
    zm1 = np.array([1.0])
    for _ in range(k - 1):
        zm1 = np.convolve(zm1, np.array([-1.0, 1.0]))
    H = np.convolve(zm1, G)[:N]
    Vc = Vcoef(U, W, k, lin)
    Vz = shift(np.pad(Vc, (0, max(0, N - len(Vc)))), -1.0)[:N]
    if len(Vz) < N:
        Vz = np.pad(Vz, (0, N - len(Vz)))
    return H - Vz, None


def unpack(q, m):
    r = np.exp(q[:m])
    return 0.5 + r * np.tanh(q[m:]), (r + 0.5) ** 2 / 8


def resid(q, k, lin, J):
    m = k - lin
    u, v = unpack(q, m)
    c0 = 2.0 ** (1 - k)
    a, _ = taylor_resid(u, v, k, lin, J + 1)
    b, _ = taylor_resid(1 - u, v, k, lin, J + 1)
    return np.concatenate([a, b]) / c0


def main():
    """continuation up the hierarchy: level-1 solution (m9_feas residual) -> level 2 -> ... -> level Jmax"""
    k, lin, Jmax, ns, seed = map(int, sys.argv[1:6])
    m = k - lin
    rng = np.random.default_rng(seed)
    if len(sys.argv) > 6:  # consistency check against the boundary equations of m9_feas
        for _ in range(3):
            u = 0.5 + rng.normal(0, .3, m)
            v = (np.abs(u - .5) + .5) ** 2 / 8 + rng.uniform(0.01, 1, m)
            a, _ = taylor_resid(u, v, k, lin, 2)
            t0, t1 = Tvals(u, v, k, lin)
            c0 = 2.0 ** (1 - k)
            print("r0 =", a[0], " c0 - V(-1) =", c0 - t0, ";  r1 =", a[1], " (1-k)c0 - V'(-1) =", (1 - k) * c0 + (k - 1) * t0 + t1)
        return
    from m9_feas import resid as resid1
    reached = []
    for it in range(ns):
        q = np.concatenate([rng.normal(0, 1.5, m), rng.normal(0, 1.0, m)])
        try:
            s = least_squares(resid1, q, args=(k, lin), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        except Exception:
            continue
        if np.abs(s.fun).max() > 1e-8:
            print(f"start {it}: level 1 not reached ({np.abs(s.fun).max():.1e})", flush=True)
            continue
        q = s.x
        top = 1
        for J in range(2, Jmax + 1):
            try:
                s = least_squares(resid, q, args=(k, lin, J), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
            except Exception:
                break
            f = np.abs(s.fun).max()
            u, v = unpack(s.x, m)
            mg = (np.sqrt(8 * v) - .5 - np.abs(u - .5)).min()
            dm = min(np.hypot(u[i] - u[j], v[i] - v[j]) for i in range(m) for j in range(i + 1, m))
            ok = f < 1e-7 and mg > 1e-6 and dm > 1e-4
            print(f"start {it}: level {J}: maxres {f:.2e} margin {mg:.2e} dmin {dm:.2e} {'OK' if ok else 'FAIL'}", flush=True)
            if not ok:
                break
            top = J
            q = s.x
            if J >= 4 or J == Jmax:
                np.set_printoptions(precision=15, linewidth=250)
                print("   u =", repr(u))
                print("   v =", repr(v))
        reached.append(top)
    print("max level reached per start:", reached)


if __name__ == "__main__":
    main()
