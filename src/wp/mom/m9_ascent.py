#!/usr/bin/env python3
"""m9 line: projected-gradient ASCENT of the order-2 quantities on the boundary set (level 1) inside the domain.

Level 1 = {four boundary equations} (resid1 of m9_feas), a manifold of dimension 2m - 4 near its regular points.
Objective: obj = min(2^k T''(-1), 2^k T~''(-1))  (necessary for level 2: both > 0, since e_2(s), e_2(alpha) > 0),
maximised by smooth surrogate  -tau log(exp(-a/tau) + exp(-b/tau))  with tangent-space projection and Gauss-Newton
retraction onto level 1.  Exact domain parametrisation (open domain <-> all of parameter space).
Reports the trajectory; a run that ends with obj > 0 exhibits a point satisfying the necessary condition.
Usage: m9_ascent.py k lin nstarts seed [steps] [min|sum] [startfile.gp]
  sum: maximise T2+ = (2^k T'' + 2^k T~'')/2 instead (e_2(s) + e_2(alpha) = 2^(k+1) T2+ must be > 0).
  startfile: gp solution file (CASE, U0, V0); starts are random perturbations of it, retracted onto level 1.
"""
import sys
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, '.')
from m9_feas import resid as resid1
from m9_ord2 import quantities


def jac(f, q, h=1e-7):
    f0 = f(q)
    J = np.zeros((len(f0), len(q)))
    for i in range(len(q)):
        e = np.zeros(len(q)); e[i] = h
        J[:, i] = (f(q + e) - f(q - e)) / (2 * h)
    return J


def retract(q, k, lin, it=30):
    for _ in range(it):
        r = resid1(q, k, lin)
        if np.abs(r).max() < 1e-11:
            return q, True
        J = jac(lambda z: resid1(z, k, lin), q)
        dq = np.linalg.lstsq(J, -r, rcond=None)[0]
        q = q + dq
        if not np.all(np.isfinite(q)):
            return q, False
    return q, np.abs(resid1(q, k, lin)).max() < 1e-9


def main():
    k, lin, ns, seed = map(int, sys.argv[1:5])
    steps = int(sys.argv[5]) if len(sys.argv) > 5 else 400
    mode = sys.argv[6] if len(sys.argv) > 6 else 'min'
    start = None
    if len(sys.argv) > 7:
        txt = open(sys.argv[7]).read()
        def grab(name):
            i = txt.index(name + ' = [')
            return np.array([float(x) for x in txt[i + len(name) + 4:txt.index(']', i)].split(',')])
        U0, V0 = grab('U0'), grab('V0')
        r0 = np.sqrt(8 * V0) - 0.5
        start = np.concatenate([np.log(r0), np.arctanh(np.clip((U0 - 0.5) / r0, -0.999999, 0.999999))])
    m = k - lin
    rng = np.random.default_rng(seed)

    def obj(q, tau=200.0):
        _, _, a, b, _, _ = quantities(q, k, lin)
        if mode == 'sum':
            return (a + b) / 2
        mn = min(a, b)
        return mn - tau * np.log(np.exp(-(a - mn) / tau) + np.exp(-(b - mn) / tau))

    for it in range(ns):
        if start is None:
            q0 = np.concatenate([rng.normal(0, 1.5, m), rng.normal(0, 1.0, m)])
            try:
                s = least_squares(resid1, q0, args=(k, lin), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
            except Exception:
                continue
            if np.abs(s.fun).max() > 1e-8:
                continue
            q = s.x
        else:
            q, ok = retract(start + (0 if it == 0 else rng.normal(0, 0.15, 2 * m)), k, lin)
            if not ok:
                continue
        eta = 0.05
        f = obj(q)
        hist = [f]
        for st in range(steps):
            J = jac(lambda z: resid1(z, k, lin), q)
            g = jac(lambda z: np.array([obj(z)]), q)[0]
            # tangent projection
            JJt = J @ J.T
            gt = g - J.T @ np.linalg.solve(JJt, J @ g)
            nrm = np.linalg.norm(gt)
            if nrm < 1e-12:
                break
            moved = False
            while eta > 1e-6:
                qn, ok = retract(q + eta * gt / nrm, k, lin)
                if ok:
                    fn = obj(qn)
                    if fn > f:
                        q, f = qn, fn
                        eta = min(eta * 1.5, 1.0)
                        moved = True
                        break
                eta /= 2
            if not moved:
                break
            hist.append(f)
        u, v, a, b, es, ea = quantities(q, k, lin)
        mg = (np.sqrt(8 * v) - .5 - np.abs(u - .5)).min()
        dm = min(np.hypot(u[i] - u[j], v[i] - v[j]) for i in range(m) for j in range(i + 1, m))
        print(f"start {it}: steps {len(hist) - 1}: obj {hist[0]:.4e} -> {f:.4e}; 2^kT'' {a:.4e} 2^kT~'' {b:.4e}; "
              f"e2(s) {es:.3e} e2(al) {ea:.3e}; margin {mg:.2e} dmin {dm:.2e} eqres {np.abs(resid1(q, k, lin)).max():.1e}", flush=True)
        np.set_printoptions(precision=6, linewidth=250)
        print("   u =", u)
        print("   v =", v)


if __name__ == "__main__":
    main()
