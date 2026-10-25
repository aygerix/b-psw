#!/usr/bin/env python3
"""ns line: independent checks of ns_core residuals.
(1) exact rational check: remainder of F - 1 mod monic g_i and of F + 1 mod monic h_i by Fraction polynomial division
    (F built by multiplying out), versus ns_core's (P0 - 1, P1) (un-normalised) at the same rational point;
    "+" side: F + 1 == 1 - P_i mod h_i, so rem(F + 1) = -(P0 - 1) - P1 t.
(2) reflection check: "+" residual of the family == "-" residual of the reflected family (a, c) -> (a, a - c), up to
    the substitution t -> -1 - t, which changes the remainder basis but not its vanishing; we compare |r| at roots.
(3) the notebook's recorded REAL k = 5 "-" half solution (b = a, d = c - a) must have small "-" residual.
"""
import sys
from fractions import Fraction as Fr
import random
import numpy as np
sys.path.insert(0, '.')
import ns_core as nc


def pmul(p, q):
    r = [Fr(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            r[i + j] += x * y
    return r


def prem(p, mon):  # coefficient lists low -> high; mon monic
    p = p[:]
    d = len(mon) - 1
    while len(p) - 1 >= d:
        lc = p[-1]
        sh = len(p) - 1 - d
        for i in range(len(mon)):
            p[sh + i] -= lc * mon[i]
        p.pop()
    return p + [Fr(0)] * (d - len(p))


def exact_rem(A, C, lin):
    fs = [[Fr(1), Fr(c) + 2, Fr(a) + Fr(c), Fr(a)] for a, c in zip(A, C)]
    if lin:
        fs.append([Fr(1), Fr(2)])
    F = [Fr(1)]
    for f in fs:
        F = pmul(F, f)
    Fm = F[:]; Fm[0] -= 1
    Fp = F[:]; Fp[0] += 1
    out = []
    for a, c in zip(A, C):
        a = Fr(a); c = Fr(c)
        g = [(c + 2) / a, (a + c) / a, Fr(1)]
        h = [Fr(2) / a, c / a, Fr(1)]
        out.append((prem(Fm, g), prem(Fp, h)))
    return out


def main():
    rng = random.Random(7)
    worst = 0.0
    for trial in range(6):
        m = rng.choice([4, 5, 6, 7, 9]); lin = rng.choice([0, 1])
        A = [Fr(rng.randint(1, 60), rng.randint(1, 20)) for _ in range(m)]
        C = [Fr(rng.randint(-60, 60), rng.randint(1, 20)) for _ in range(m)]
        ex = exact_rem(A, C, lin)
        a = np.array([[float(x) for x in A]], dtype=complex); c = np.array([[float(x) for x in C]], dtype=complex)
        for side, k in (('-', 0), ('+', 1)):
            if side == '-':
                p = (a + c) / a; q = (c + 2) / a
            else:
                p = c / a; q = 2 / a
            x0, x1 = nc.reduce_members(a, c, p, q)
            idx = np.arange(m); x0[:, idx, idx] = 1; x1[:, idx, idx] = 0
            P0 = np.ones((1, m), complex); P1 = np.zeros((1, m), complex)
            for j in range(m):
                y0 = x0[:, :, j]; y1 = x1[:, :, j]
                P0, P1 = P0 * y0 - q * P1 * y1, P0 * y1 + P1 * y0 - p * P1 * y1
            if lin:
                P0, P1 = P0 - q * P1 * 2, P0 * 2 + P1 - p * P1 * 2
            for i in range(m):
                e0, e1 = [float(z) for z in ex[i][k]]
                n0, n1 = (P0[0, i] - 1).real, P1[0, i].real
                if side == '+':
                    n0, n1 = -n0, -n1
                rel = max(abs(e0 - n0), abs(e1 - n1)) / (1 + abs(e0) + abs(e1))
                worst = max(worst, rel)
        print(f"trial {trial}: m = {m} lin = {lin}: worst relative discrepancy so far {worst:.2e}")
    # (2) "+" side by direct evaluation at the roots of h_i, and versus the "-" side of the reflected family
    rs = np.random.default_rng(3)
    pass
    for lin in (0, 1):
        mm = 7 - lin   # k = mm + lin = 7 odd (the sign identity P~(w) = (-1)^(k-1) P(z) needs k odd)
        a = rs.uniform(0.5, 5, (4, mm)) + 0j; c = rs.uniform(-1, 1, (4, mm)) * np.sqrt(8 * a.real) * 0.5 + 0j
        e0, e1, p, q = nc.side_resid(a, c, lin, '+', raw=True)
        z = (-p + np.sqrt(p * p - 4 * q + 0j)) / 2           # a root of h_i
        val_rem = e0 + e1 * z
        F = lambda t, aa, cc: 2 * t + 1 + t * (t + 1) * (aa * t + cc)
        val_dir = np.ones_like(z)
        for j in range(mm):
            fj = F(z, a[:, j:j+1], c[:, j:j+1])
            fj[:, j] = 1.0
            val_dir = val_dir * fj
        if lin:
            val_dir = val_dir * (2 * z + 1)
        val_dir = val_dir - 1
        f0, f1, pr, qr = nc.side_resid(a, a - c, lin, '-', raw=True)
        w = -1 - z                                            # root of the reflected g_i
        val_ref = f0 + f1 * w
        print(f"lin={lin}: '+' remainder vs direct evaluation: max rel err",
              (np.abs(val_rem - val_dir) / (1 + np.abs(val_dir))).max(),
              "; vs reflected '-' : max rel err", (np.abs(val_ref - val_dir) / (1 + np.abs(val_dir))).max())
    # (3) notebook k = 5 real "-" half solution
    b = np.array([0.517715, 1.990854, 21.097857, 0.896437, 6.526845])
    d = np.array([-1.705746, -3.59832, -6.677547, -2.316047, -6.113322])
    a = (b + 0j)[None]; c = (b + d + 0j)[None]
    r0, r1 = nc.side_resid(a, c, 0, '-')
    s0, s1 = nc.side_resid(a, c, 0, '+')
    print("notebook k=5 '-' half solution: '-' residual max", np.abs(np.concatenate([r0, r1], 1)).max(),
          " '+' residual max", np.abs(np.concatenate([s0, s1], 1)).max())


if __name__ == "__main__":
    main()
