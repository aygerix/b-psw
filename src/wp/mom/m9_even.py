#!/usr/bin/env python3
"""m9 line: relaxed EVEN moment problem behind the order-2 obstruction at k = 9.

Exact facts (gp: m9_mform.out, m9_t2.out), with S = M1_0, q = M1_1^2, moments over the m cubic members:
  9 cubics:   Z0+ = 39/8 S^2 + S(M12 - 2M20 - 321/32) + q/2 - 51/4 M12 + 261/8 M20 - 9/2 M14 + 57/2 M22 - 20 M30 + 255/256
              Z1+ = -73/2 S^2 + S(-8M12 + 16M20 + 933/16) - 4q + 90 M12 - 234 M20 + 35 M14 - 222 M22 + 156 M30
              T2+ = Psi := S^2 - 209/4 S - 11 M12 + 23 M20 + 1785/32         on {Z0+ = Z1+ = 0}
  8 + 2t+1:   Z0+ = 31/8 S^2 + S(M12 - 2M20 - 6) + q/2 - 9 M12 + 191/8 M20 - 4 M14 + 51/2 M22 - 18 M30 + 127/256
              Z1+ = -115/4 S^2 + S(-8M12 + 16M20 + 33) - 4q + 249/4 M12 - 675/4 M20 + 31 M14 - 198 M22 + 140 M30 + 1/2
              Psi := S^2 - 40 S - 10 M12 + 21 M20 + 1113/32
A real cubic family needs Psi > 0 (Psi = (e2(s) + e2(alpha)) / 2^(k+1), s, alpha > 0).
RELAXATION here: atoms (v_i, y_i = w_i^2) in the closed domain v >= 1/32, 0 <= y <= (sqrt(8v) - 1/2)^2; the odd
equations are dropped; q is eliminated (8 Z0+ + Z1+ is q-free) and then only 0 <= q <= S M12 (Cauchy-Schwarz) is kept
(mode 'q') or dropped (mode 'noq').  Maximise Psi by projected gradient ascent on the constraint manifold.
Usage: m9_even.py lin mode nstarts seed
"""
import sys
import numpy as np


def moments(p, m):
    r = np.exp(np.clip(p[:m], -40, 4))
    v = (r + .5) ** 2 / 8
    y = (r * (1 + np.tanh(p[m:2 * m])) / 2) ** 2
    M = {'S': v.sum(), 'M12': (v * y).sum(), 'M20': (v * v).sum(), 'M14': (v * y * y).sum(),
         'M22': (v * v * y).sum(), 'M30': (v ** 3).sum()}
    return M, v, y


def forms(M, lin):
    S, M12, M20, M14, M22, M30 = M['S'], M['M12'], M['M20'], M['M14'], M['M22'], M['M30']
    if lin == 0:
        z0 = 39 / 8 * S ** 2 + S * (M12 - 2 * M20 - 321 / 32) - 51 / 4 * M12 + 261 / 8 * M20 - 4.5 * M14 + 28.5 * M22 - 20 * M30 + 255 / 256
        z1 = -36.5 * S ** 2 + S * (-8 * M12 + 16 * M20 + 933 / 16) + 90 * M12 - 234 * M20 + 35 * M14 - 222 * M22 + 156 * M30
        psi = S ** 2 - 209 / 4 * S - 11 * M12 + 23 * M20 + 1785 / 32
    else:
        z0 = 31 / 8 * S ** 2 + S * (M12 - 2 * M20 - 6) - 9 * M12 + 191 / 8 * M20 - 4 * M14 + 25.5 * M22 - 18 * M30 + 127 / 256
        z1 = -115 / 4 * S ** 2 + S * (-8 * M12 + 16 * M20 + 33) + 249 / 4 * M12 - 675 / 4 * M20 + 31 * M14 - 198 * M22 + 140 * M30 + 0.5
        psi = S ** 2 - 40 * S - 10 * M12 + 21 * M20 + 1113 / 32
    # z0 + q/2 = 0 and z1 - 4q = 0  ->  q = -2 z0 ; q-free equation 8 z0 + z1 = 0
    return 8 * z0 + z1, -2 * z0, psi


def main():
    """variables: atoms (rho_i, z_i) and zeta; q = theta * S * M12, theta = (1 + tanh zeta)/2 in (0, 1)  [mode 'q'];
    mode 'noq': q free real (the unbounded relaxation)."""
    lin, mode, ns, seed = int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    m = 9 - lin
    rng = np.random.default_rng(seed)

    def parts(p):
        M, v, y = moments(p[:2 * m], m)
        S, M12, M20, M14, M22, M30 = M['S'], M['M12'], M['M20'], M['M14'], M['M22'], M['M30']
        q = (1 + np.tanh(p[2 * m])) / 2 * S * M12 if mode == 'q' else p[2 * m]
        if lin == 0:
            z0 = 39 / 8 * S ** 2 + S * (M12 - 2 * M20 - 321 / 32) + q / 2 - 51 / 4 * M12 + 261 / 8 * M20 - 4.5 * M14 + 28.5 * M22 - 20 * M30 + 255 / 256
            z1 = -36.5 * S ** 2 + S * (-8 * M12 + 16 * M20 + 933 / 16) - 4 * q + 90 * M12 - 234 * M20 + 35 * M14 - 222 * M22 + 156 * M30
            psi = S ** 2 - 209 / 4 * S - 11 * M12 + 23 * M20 + 1785 / 32
        else:
            z0 = 31 / 8 * S ** 2 + S * (M12 - 2 * M20 - 6) + q / 2 - 9 * M12 + 191 / 8 * M20 - 4 * M14 + 25.5 * M22 - 18 * M30 + 127 / 256
            z1 = -115 / 4 * S ** 2 + S * (-8 * M12 + 16 * M20 + 33) - 4 * q + 249 / 4 * M12 - 675 / 4 * M20 + 31 * M14 - 198 * M22 + 140 * M30 + 0.5
            psi = S ** 2 - 40 * S - 10 * M12 + 21 * M20 + 1113 / 32
        return np.array([z0, z1]), psi, M, v, y, q

    def jac(f, p, h=1e-7):
        f0 = np.atleast_1d(f(p))
        J = np.zeros((len(f0), len(p)))
        for i in range(len(p)):
            e = np.zeros(len(p)); e[i] = h
            J[:, i] = (np.atleast_1d(f(p + e)) - np.atleast_1d(f(p - e))) / (2 * h)
        return J

    cons = lambda p: parts(p)[0]
    obj = lambda p: parts(p)[1]

    def retract(p):
        for _ in range(60):
            c = cons(p)
            if np.abs(c).max() < 1e-11:
                return p, True
            J = jac(cons, p)
            if not np.all(np.isfinite(J)) or not np.all(np.isfinite(c)):
                return p, False
            try:
                p = p - np.linalg.lstsq(J, c, rcond=None)[0]
            except np.linalg.LinAlgError:
                return p, False
            if not np.all(np.isfinite(p)):
                return p, False
        return p, np.abs(cons(p)).max() < 1e-9

    best = None
    for it in range(ns):
        p = np.concatenate([rng.normal(-0.5, 1.2, m), rng.normal(0, 1.5, m), [rng.normal(0, 1)]])
        p, ok = retract(p)
        if not ok:
            continue
        f = obj(p)
        eta = 0.1
        for st in range(4000):
            g = jac(obj, p)[0]
            J = jac(cons, p)
            try:
                gt = g - J.T @ np.linalg.solve(J @ J.T, J @ g)
            except np.linalg.LinAlgError:
                break
            nr = np.linalg.norm(gt)
            if nr < 1e-10:
                break
            moved = False
            while eta > 1e-9:
                pn, ok = retract(p + eta * gt / nr)
                if ok and obj(pn) > f:
                    p, f = pn, obj(pn)
                    eta *= 1.5
                    moved = True
                    break
                eta /= 2
            if not moved:
                break
        c, psi, M, v, y, q = parts(p)
        np.set_printoptions(precision=5, linewidth=250)
        print(f"start {it}: Psi {psi:.6f}  (eqs {np.abs(c).max():.1e}, q {q:.5f}, S*M12 {M['S'] * M['M12']:.5f}, S {M['S']:.5f}, M20 {M['M20']:.5f}, M12 {M['M12']:.5f})", flush=True)
        print("   v =", v, "\n   y/ymax =", y / (np.sqrt(8 * v) - .5) ** 2)
        if best is None or psi > best:
            best = psi
    print("sup Psi found:", best)


if __name__ == "__main__":
    main()
