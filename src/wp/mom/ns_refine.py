#!/usr/bin/env python3
"""ns line: high-precision refinement / verification of a candidate cubic Williams family (mpmath).
Residual: coefficients of  (F - 1) mod monic(t G)  and  (F + 1) mod monic((t+1) H)  (G = prod g_i, H = prod h_i),
an exact reformulation of the family condition when the g_i (resp. h_i) are pairwise coprime and prime to t (t+1).
Gauss-Newton (least squares, QR via mpmath) in the 2m (real) or 2m complex unknowns (a_i, c_i).
Usage: python3 ns_refine.py <npz file> <index> [dps=60] [iters=60] [lin] [sides=-+]
       or import refine(a, c, lin, ...) .
Prints the residual per iteration, the smallest singular value of the Jacobian, member separation and margins.
"""
import sys
import mpmath as mp


def pmul(p, q):
    r = [mp.mpc(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        if x == 0:
            continue
        for j, y in enumerate(q):
            r[i + j] += x * y
    return r


def prem(p, mon):
    """remainder of p modulo monic polynomial mon (coefficient lists, low -> high)"""
    p = list(p); d = len(mon) - 1
    while len(p) - 1 >= d:
        lc = p[-1]; sh = len(p) - 1 - d
        if lc != 0:
            for i in range(d):
                p[sh + i] -= lc * mon[i]
        p.pop()
    return p + [mp.mpc(0)] * (d - len(p))


def residual(a, c, lin, sides='-+'):
    m = len(a)
    fs = [[mp.mpc(1), c[i] + 2, a[i] + c[i], a[i]] for i in range(m)]
    if lin:
        fs.append([mp.mpc(1), mp.mpc(2)])
    F = [mp.mpc(1)]
    for f in fs:
        F = pmul(F, f)
    out = []
    for s in sides:
        if s == '-':
            mods = [[(c[i] + 2) / a[i], (a[i] + c[i]) / a[i], mp.mpc(1)] for i in range(m)]
            M = [mp.mpc(0), mp.mpc(1)]          # t
            T = list(F); T[0] -= 1
        else:
            mods = [[2 / a[i], c[i] / a[i], mp.mpc(1)] for i in range(m)]
            M = [mp.mpc(1), mp.mpc(1)]          # t + 1
            T = list(F); T[0] += 1
        for q in mods:
            M = pmul(M, q)
        r = prem(T, M)
        # scale: divide by the size of F's coefficients so the residual is relative
        sc = max(abs(x) for x in F)
        out += [x / sc for x in r]
    return out


def refine(a, c, lin, sides='-+', real=True, dps=60, iters=60, verbose=True):
    mp.mp.dps = dps
    m = len(a)
    if real:
        x = [mp.mpf(z.real if isinstance(z, complex) else z) for z in list(a) + list(c)]
    else:
        x = [mp.mpc(z) for z in list(a) + list(c)]
    hist = []
    for it in range(iters):
        A = [mp.mpc(z) for z in x[:m]]; C = [mp.mpc(z) for z in x[m:]]
        r = residual(A, C, lin, sides)
        nr = max(abs(z) for z in r)
        hist.append(nr)
        if verbose:
            print(f"   it {it}: max |residual| = {mp.nstr(nr, 5)}", flush=True)
        if nr < mp.mpf(10) ** (-(dps - 8)):
            break
        # Jacobian by complex-step-free finite differences at high precision
        h = mp.mpf(10) ** (-(dps // 2))
        J = []
        for j in range(2 * m):
            x2 = list(x); x2[j] += h
            A2 = [mp.mpc(z) for z in x2[:m]]; C2 = [mp.mpc(z) for z in x2[m:]]
            r2 = residual(A2, C2, lin, sides)
            J.append([(u - v) / h for u, v in zip(r2, r)])
        nrow = len(r)
        if real:
            Jm = mp.matrix(2 * nrow, 2 * m); rv = mp.matrix(2 * nrow, 1)
            for i in range(nrow):
                rv[i] = r[i].real; rv[nrow + i] = r[i].imag
                for j in range(2 * m):
                    Jm[i, j] = J[j][i].real; Jm[nrow + i, j] = J[j][i].imag
        else:
            Jm = mp.matrix(nrow, 2 * m); rv = mp.matrix(nrow, 1)
            for i in range(nrow):
                rv[i] = r[i]
                for j in range(2 * m):
                    Jm[i, j] = J[j][i]
        JH = Jm.H if not real else Jm.T
        dx = mp.lu_solve(JH * Jm, JH * rv)
        x = [x[j] - dx[j] for j in range(2 * m)]
    sv = mp.svd(Jm, compute_uv=False)
    smin = min(sv[i] for i in range(len(sv)))
    A = [mp.mpc(z) for z in x[:m]]; C = [mp.mpc(z) for z in x[m:]]
    return A, C, hist, smin


def main():
    import numpy as np
    f, i = sys.argv[1], int(sys.argv[2])
    dps = int(sys.argv[3]) if len(sys.argv) > 3 else 60
    iters = int(sys.argv[4]) if len(sys.argv) > 4 else 60
    lin = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    sides = sys.argv[6] if len(sys.argv) > 6 else '-+'
    d = np.load(f)
    a, c = d['a'][i], d['c'][i]
    real = np.abs(a.imag).max() + np.abs(c.imag).max() == 0
    A, C, hist, smin = refine(list(a), list(c), lin, sides, real, dps, iters)
    print("final residual", mp.nstr(hist[-1], 5), " Jacobian smallest singular value", mp.nstr(smin, 5))
    print("a =", [mp.nstr(z, 20) for z in A]); print("c =", [mp.nstr(z, 20) for z in C])


if __name__ == "__main__":
    main()
