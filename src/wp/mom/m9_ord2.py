#!/usr/bin/env python3
"""m9 line: the order-2 necessary condition on top of the four boundary equations.

On the boundary set (T(-1) = 2^(1-k), T'(-1) = 0) the order-2 Taylor equation at t = -1 reads
    e_2(s) = 2^k T''(-1),     s_i = f_i'(-1) = h_i(-1) = 2 + (1 - u_i)/v_i > 0   (h_i > 0: no real root, a_i > 0),
(derivation: V''(-1) = (k-1)(k-2) T(-1) + (2k-2) T'(-1) + T''(-1) and e_2(s) = 2^k V''(-1) - 2(k-1)(k-2)),
and the same for the reflected family (alpha_i = 2 + u_i/v_i).  Linear member: s = alpha = 2.
NECESSARY (since e_2 > 0): T''(-1) > 0 and T~''(-1) > 0.
Modes:
  maxmin : maximise min(2^k T''(-1), 2^k T~''(-1)) over the boundary set inside the domain (SLSQP, exact domain
           parametrisation, 4 equality constraints).  A negative supremum would mean level 2 is infeasible.
  level2 : solve the 6 equations (boundary + both order-2 equations) by least squares from random starts.
Usage: m9_ord2.py k lin mode nstarts seed
"""
import sys
import numpy as np
from scipy.optimize import least_squares, minimize
sys.path.insert(0, '.')
from m9_feas import resid as resid1


def Tcoef(U, W, k, lin):
    T = np.zeros(k)
    T[0] = 1.0
    for u, v in zip(U, W):
        d = np.zeros(k)
        d[0] = 1.0
        d[1] = -(1 + u)
        for j in range(2, k):
            d[j] = -(1 + u) * d[j - 1] - (u + 2 * v) * d[j - 2]
        f = np.zeros(k)
        f[0] = 1.0
        f[3:] = v * d[:k - 3]
        T = np.convolve(T, f)[:k]
    if lin:
        T = np.convolve(T, np.array([1.0, 0.5] + [0.0] * (k - 2)))[:k]
    return T


def T2(U, W, k, lin):
    T = Tcoef(U, W, k, lin)
    j = np.arange(k)
    return np.dot(j * (j - 1) * T, (-1.0) ** j)


def unpack(q, m):
    r = np.exp(q[:m])
    return 0.5 + r * np.tanh(q[m:]), (r + 0.5) ** 2 / 8


def e2(x):
    return (x.sum() ** 2 - (x ** 2).sum()) / 2


def quantities(q, k, lin):
    m = k - lin
    u, v = unpack(q, m)
    s, al = 2 + (1 - u) / v, 2 + u / v
    if lin:
        s, al = np.append(s, 2.0), np.append(al, 2.0)
    return u, v, 2.0 ** k * T2(u, v, k, lin), 2.0 ** k * T2(1 - u, v, k, lin), e2(s), e2(al)


def main():
    k, lin, mode, ns, seed = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
    m = k - lin
    rng = np.random.default_rng(seed)
    best = None
    for it in range(ns):
        q0 = np.concatenate([rng.normal(0, 1.5, m), rng.normal(0, 1.0, m)])
        try:
            s = least_squares(resid1, q0, args=(k, lin), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        except Exception:
            continue
        if np.abs(s.fun).max() > 1e-8:
            continue
        q = s.x
        if mode == 'maxmin':
            # maximise tau s.t. boundary eqs = 0, 2^k T'' >= tau, 2^k T~'' >= tau
            def cin(z):
                _, _, a, b, _, _ = quantities(z[:-1], k, lin)
                return np.array([a - z[-1], b - z[-1]]) / 1e3
            _, _, a0, b0, _, _ = quantities(q, k, lin)
            z0 = np.append(q, min(a0, b0))
            try:
                res = minimize(lambda z: -z[-1] / 1e3, z0, method='SLSQP',
                               constraints=[{'type': 'eq', 'fun': lambda z: resid1(z[:-1], k, lin)},
                                            {'type': 'ineq', 'fun': cin}],
                               options={'maxiter': 2000, 'ftol': 1e-12})
                z = res.x
            except Exception as e:
                print("exception", e)
                z = z0
            er = np.abs(resid1(z[:-1], k, lin)).max()
            u, v, a, b, es, ea = quantities(z[:-1], k, lin)
            mg = (np.sqrt(8 * v) - .5 - np.abs(u - .5)).min()
            val = min(a, b)
            print(f"start {it}: start min(2^k T'', 2^k T~'') = {min(a0, b0):.4e} -> {val:.4e}  (eqres {er:.1e}, margin {mg:.2e}); "
                  f"e2(s) {es:.3e}, e2(alpha) {ea:.3e}", flush=True)
            if er < 1e-7 and mg > 0 and (best is None or val > best[0]):
                best = (val, u, v, a, b, es, ea)
        elif mode.startswith('theta'):
            th = float(mode[5:])
            def r6t(z):
                u, v, a, b, es, ea = quantities(z, k, lin)
                sc = 2.0 ** (k - 1)
                return np.concatenate([resid1(z, k, lin), [(a - th) / sc, (b - th) / sc]])
            try:
                s2 = least_squares(r6t, q, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=6000)
            except Exception:
                continue
            u, v, a, b, es, ea = quantities(s2.x, k, lin)
            mg = (np.sqrt(8 * v) - .5 - np.abs(u - .5)).min()
            dm = min(np.hypot(u[i] - u[j], v[i] - v[j]) for i in range(m) for j in range(i + 1, m))
            f = np.abs(s2.fun).max()
            ok = f < 1e-7 and mg > 1e-6 and dm > 1e-4
            print(f"start {it}: theta {th}: maxres {f:.2e} margin {mg:.2e} dmin {dm:.2e} {'OK' if ok else 'FAIL'}; 2^kT'' {a:.4e} 2^kT~'' {b:.4e}; e2(s) {es:.3e} e2(al) {ea:.3e}", flush=True)
            if ok:
                np.set_printoptions(precision=15, linewidth=250)
                print("   u =", repr(u))
                print("   v =", repr(v))
        else:
            def r6(z):
                u, v, a, b, es, ea = quantities(z, k, lin)
                sc = 2.0 ** (k - 1)
                return np.concatenate([resid1(z, k, lin), [(a - es) / sc, (b - ea) / sc]])
            try:
                s2 = least_squares(r6, q, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=6000)
            except Exception:
                continue
            u, v, a, b, es, ea = quantities(s2.x, k, lin)
            mg = (np.sqrt(8 * v) - .5 - np.abs(u - .5)).min()
            f = np.abs(s2.fun).max()
            print(f"start {it}: level-2 maxres {f:.2e} margin {mg:.2e}; 2^kT'' {a:.3e} e2(s) {es:.3e}; 2^kT~'' {b:.3e} e2(al) {ea:.3e}", flush=True)
    if best is not None:
        np.set_printoptions(precision=15, linewidth=250)
        print("BEST min(2^k T'', 2^k T~'') =", best[0], " values", best[3], best[4], " e2(s), e2(alpha):", best[5], best[6])
        print("u =", repr(best[1]))
        print("v =", repr(best[2]))


if __name__ == "__main__":
    main()
