#!/usr/bin/env python3
"""m9 line: does the full level-1 set (four boundary equations) reach T2+ > 0 at LARGE S (members with large v)?
Start from the even-relaxation point with Psi > 0 (m9_even.py 1 q, start 13), attach random signs w = +-sqrt(y),
solve the 4 boundary equations by least squares nearby, then report 2^k T'', 2^k T~'' and e2(s), e2(alpha).
Usage: m9_bigS.py ntrials seed
"""
import sys
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, '.')
from m9_feas import resid as resid1
from m9_ord2 import quantities
k, lin = 9, 1
m = 8
v0 = np.array([5.31592e+01, 3.66598e-02, 1.72504e+00, 6.03071e-01, 2.31615e-01, 1.45983e+02, 4.64048e-02, 2.05488e-01])
fr = np.array([0.95518, 0.7409, 0.0587, 0.76496, 0.04993, 0.11507, 0.01456, 0.92046])
nt, seed = map(int, sys.argv[1:3])
rng = np.random.default_rng(seed)
for it in range(nt):
    r = np.sqrt(8 * v0) - .5
    w = r * np.sqrt(fr) * rng.choice([-1, 1], m)
    q0 = np.concatenate([np.log(r), np.arctanh(np.clip(w / r, -.9999, .9999))]) + rng.normal(0, 0.05, 2 * m)
    s = least_squares(resid1, q0, args=(k, lin), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=5000)
    f = np.abs(s.fun).max()
    u, v, a, b, es, ea = quantities(s.x, k, lin)
    mg = (np.sqrt(8 * v) - .5 - np.abs(u - .5)).min()
    print(f"trial {it}: level1 res {f:.1e}  S {v.sum():.3f} vmax {v.max():.3f}  2^kT'' {a:.4e} 2^kT~'' {b:.4e} e2(s) {es:.3e} e2(al) {ea:.3e} margin {mg:.2e}", flush=True)
