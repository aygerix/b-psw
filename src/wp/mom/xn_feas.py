#!/usr/bin/env python3
"""xn line: numerical feasibility of the order-N truncation (family + reflection: 2N equations) on the real domain.

Domain (exact parametrisation): r = exp(rho) > 0, v = (r + 1/2)^2 / 8, u = 1/2 + r tanh(z)  <=>  (|u-1/2|+1/2)^2 < 8v.
Equations: E_n = 0 (n < N) for the family and the reflected family (xn_core.E_both), scaled by 2^(k-1).
N = 2: the four boundary equations; N = 2k: the full Williams condition.
Usage: xn_feas.py k lin N nstarts seed [maxnfev]
COMPUTED evidence only: failure to converge is not a proof of infeasibility.
"""
import sys
import numpy as np
from scipy.optimize import least_squares
from xn_core import E_both


def unpack(p):
    m = len(p) // 2
    r = np.exp(p[:m])
    v = (r + 0.5) ** 2 / 8
    u = 0.5 + r * np.tanh(p[m:])
    return u, v


def resid(p, k, lin, N):
    u, v = unpack(p)
    e = E_both([float(x) for x in u], [float(x) for x in v], lin, N)
    return np.array(e) * 2.0 ** (k - 1)


def margins(u, v):
    mg = np.sqrt(8 * v) - (np.abs(u - 0.5) + 0.5)
    pts = np.stack([u, v], 1)
    dmin = min(np.linalg.norm(pts[i] - pts[j]) for i in range(len(u)) for j in range(i + 1, len(u))) if len(u) > 1 else 1.0
    return mg.min(), dmin


def main():
    k, lin, N, ns, seed = map(int, sys.argv[1:6])
    mx = int(sys.argv[6]) if len(sys.argv) > 6 else 3000
    m = k - lin
    rng = np.random.default_rng(seed)
    nok, best = 0, np.inf
    for it in range(ns):
        p0 = np.concatenate([rng.normal(0, 1.5, m), rng.normal(0, 1.0, m)])
        try:
            with np.errstate(all="ignore"):
                sol = least_squares(resid, p0, args=(k, lin, N), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=mx)
        except Exception:
            continue
        if not np.all(np.isfinite(sol.fun)):
            continue
        f = np.abs(sol.fun).max()
        best = min(best, f)
        u, v = unpack(sol.x)
        mg, dm = margins(u, v)
        tag = "OK  " if f < 1e-9 else "fail"
        print(f"start {it}: {tag} maxres {f:.3e} margin {mg:.3e} dmin {dm:.3e} S {v.sum():.6f} maxv {v.max():.3e}", flush=True)
        if f < 1e-9:
            nok += 1
            print("   u =", np.array2string(u, precision=10, max_line_width=300))
            print("   v =", np.array2string(v, precision=10, max_line_width=300))
    print(f"k={k} lin={lin} N={N}: converged {nok} of {ns}; best maxres {best:.3e}")


if __name__ == "__main__":
    main()
