#!/usr/bin/env python3
"""aud13: discriminating numerical feasibility test of the xn cover on ACTUAL member configurations (k = 9).
For every certificate and sampled s0 in its S-interval, minimise the residual of the constraints that the certificate
uses (E0 = E1 = 0 always; PHI >= 0, POS(t0) >= 0 as hinges; CS holds automatically for real configurations) over
configurations with S = s0 exactly:  v = 1/32 + (s0 - m/32) softmax(alpha), w = (sqrt(8v) - 1/2) sin(theta).
A residual ~ 0 at any s0 would contradict the certificate (refute the cover or expose a checker bug).
Rows: aud13_rows.json via aud13_cert.  Usage: aud13_opt.py lin nstarts seed [first last]
"""
import sys

import numpy as np
from scipy.optimize import least_squares

import aud13_cert as AC

lin, nst, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
m = 9 - lin
rng = np.random.default_rng(seed)
CERTS = [AC.parse_gp(l.strip()) for l in open(f"xn_cert9_l{lin}.gp") if l.strip()]
first = int(sys.argv[4]) if len(sys.argv) > 4 else 1
last = int(sys.argv[5]) if len(sys.argv) > 5 else len(CERTS)


def fform(name):
    col, pc, c = AC.form(lin, name)
    return [(a, j, float(c0), float(c1)) for (a, j), (c0, c1) in col.items()], float(pc), [float(x) for x in c]


FORMS = {nm: fform(nm) for nm in ["E0", "E1", "PHI"] + [f"T{n}" for n in range(9)]}


def ev(name, v, w):
    cols, pc, c = FORMS[name]
    S = v.sum()
    W = w * w
    val = sum((c0 + c1 * S) * np.sum(v ** a * W ** j) for a, j, c0, c1 in cols)
    P = np.sum(v * w) ** 2 / 2
    return val + pc * P + c[0] + c[1] * S + c[2] * S * S


def pos(t0, v, w):
    return sum(ev(f"T{n}", v, w) * t0 ** (8 - n) for n in range(9)) - t0 ** 8 / 256


def members(x, s0):
    al, th = x[:m], x[m:]
    e = np.exp(al - al.max())
    v = 1 / 32 + (s0 - m / 32) * e / e.sum()
    w = (np.sqrt(8 * v) - 0.5) * np.sin(th)
    return v, w


def resid(x, s0, kinds):
    v, w = members(x, s0)
    r = [ev("E0", v, w), ev("E1", v, w)]
    for k, p in kinds:
        g = ev("PHI", v, w) if k == "PHI" else pos(p, v, w)
        r.append(min(g, 0.0))
    return np.array(r)


worst = []
for ci in range(first - 1, last):
    c = CERTS[ci]
    Sa, Sb = float(AC.SofR(c[1], m)), float(AC.SofR(c[2], m))
    kinds = [(k, float(p)) for k, p, z in c[4] if k in ("PHI", "POS")]
    if len(sys.argv) > 6 and sys.argv[6] == "noineq":     # positive control: drop the inequality cuts
        kinds = []
    best = None
    for it in range(nst):
        s0 = rng.uniform(Sa, Sb)
        x0 = np.concatenate([rng.normal(0, rng.uniform(0.3, 3), m), rng.uniform(-np.pi / 2, np.pi / 2, m)])
        try:
            res = least_squares(resid, x0, args=(s0, kinds), method="trf", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=int(__import__("os").environ.get("NFEV", "3000")))
        except Exception:
            continue
        nr = float(np.linalg.norm(res.fun))
        if best is None or nr < best[0]:
            best = (nr, s0, res.x)
    v, w = members(best[2], best[1])
    print(f"cert {ci + 1:2d} lin={lin} S=[{Sa:.4f},{Sb:.4f}] rows {[k for k, _ in kinds] or ['eq only']}: min residual {best[0]:.3e} "
          f"at S={best[1]:.4f} (E0, E1 = {ev('E0', v, w):.2e}, {ev('E1', v, w):.2e}; Phi = {ev('PHI', v, w):.3g})", flush=True)
    worst.append(best[0])
print(f"lin={lin}: certificates {first}..{last}: smallest min-residual over all {min(worst):.3e} "
      f"({sum(1 for x in worst if x < 1e-8)} below 1e-8 -> each would be a contradiction)")
