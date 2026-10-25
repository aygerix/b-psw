#!/usr/bin/env python3
"""m9 line: Newton-homotopy climb of the Taylor hierarchy.

From a level-(J-1) point x0 (orders 0..J-1 at t = -1 vanish for both halves), track
    G_lam(x) = [ r_0..r_(J-1) (both halves) ;  r_J(x) - (1 - lam) r_J(x0) (both halves) ] = 0,   lam: 0 -> 1,
with minimum-norm Gauss-Newton corrections (underdetermined system, 2m unknowns) in the exact domain parametrisation.
Success: lam = 1 reached with residual < tol, margin > 1e-6, member distance > 1e-4.  Then continue to J + 1.
Usage: m9_homot.py k lin J0 Jmax seed 'u,..' 'v,..' [nsteps] [tol]
"""
import sys
import numpy as np
sys.path.insert(0, '.')
from m9_hier import taylor_resid, unpack

k, lin, J0, Jmax, seed = map(int, sys.argv[1:6])
U = np.array([float(x) for x in sys.argv[6].split(',')]); V = np.array([float(x) for x in sys.argv[7].split(',')])
nsteps = int(sys.argv[8]) if len(sys.argv) > 8 else 40
TOL = float(sys.argv[9]) if len(sys.argv) > 9 else 1e-10
m = k - lin
c0 = 2.0 ** (1 - k)
r0 = np.sqrt(8 * V) - .5
q = np.concatenate([np.log(r0), np.arctanh(np.clip((U - .5) / r0, -.999999, .999999))])


def rvec(q, N):
    u, v = unpack(q, m)
    a, _ = taylor_resid(u, v, k, lin, N)
    b, _ = taylor_resid(1 - u, v, k, lin, N)
    return np.concatenate([a, b]) / c0   # [fam r_0..r_(N-1), refl r_0..r_(N-1)]


def G(q, J, target):
    r = rvec(q, J + 1)
    fam, ref = r[:J + 1], r[J + 1:]
    return np.concatenate([fam[:J], ref[:J], [fam[J] - target[0], ref[J] - target[1]]])


def jac(f, q, h=1e-7):
    f0 = f(q)
    Jm = np.zeros((len(f0), len(q)))
    for i in range(len(q)):
        e = np.zeros(len(q)); e[i] = h
        Jm[:, i] = (f(q + e) - f(q - e)) / (2 * h)
    return Jm


def correct(q, f, it=30, tol=None):
    tol = TOL if tol is None else tol
    for _ in range(it):
        r = f(q)
        if not np.all(np.isfinite(r)):
            return q, False
        if np.abs(r).max() < tol:
            return q, True
        Jm = jac(f, q)
        dq = np.linalg.lstsq(Jm, -r, rcond=None)[0]
        nd = np.linalg.norm(dq)
        if nd > 0.5:
            dq *= 0.5 / nd
        q = q + dq
    r = f(q)
    return q, np.all(np.isfinite(r)) and np.abs(r).max() < tol * 100


np.set_printoptions(precision=16, linewidth=300)
for J in range(J0, Jmax + 1):
    r = rvec(q, J + 1)
    rJ0 = np.array([r[J], r[2 * J + 1]])
    lam, dl, ok = 0.0, 1.0 / nsteps, True
    while lam < 1.0 - 1e-12:
        lt = min(1.0, lam + dl)
        tgt = (1 - lt) * rJ0
        qn, good = correct(q, lambda z: G(z, J, tgt))
        u, v = unpack(qn, m)
        mg = (np.sqrt(8 * v) - .5 - np.abs(u - .5)).min() if good else -1
        dm = min(np.hypot(u[i] - u[j], v[i] - v[j]) for i in range(m) for j in range(i + 1, m)) if good else -1
        if good and mg > 1e-6 and dm > 1e-4:
            q, lam = qn, lt
            dl = min(dl * 1.5, 0.2)
        else:
            dl /= 2
            if dl < 1e-6:
                ok = False
                break
    u, v = unpack(q, m)
    mg = (np.sqrt(8 * v) - .5 - np.abs(u - .5)).min()
    dm = min(np.hypot(u[i] - u[j], v[i] - v[j]) for i in range(m) for j in range(i + 1, m))
    print(f"level {J}: {'REACHED' if ok else 'STUCK at lam = %.6f' % lam}; margin {mg:.3e} dmin {dm:.3e} S {v.sum():.3f} vmax {v.max():.3f}"
          f" amin {1 / v.max():.4f}; residual {np.abs(G(q, J, (1 - lam) * rJ0)).max():.1e}", flush=True)
    print("   u =", repr(u)); print("   v =", repr(v), flush=True)
    if not ok:
        break
