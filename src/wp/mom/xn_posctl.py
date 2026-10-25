#!/usr/bin/env python3
"""xn line: POSITIVE CONTROL for the Farkas certificates.

Solutions of the four boundary equations (N = 2) exist numerically for k = 9.  Each one satisfies every constraint
of the relaxation except possibly PHI / N3 / POS.  Hence no such solution may lie in an S-interval certified by a
certificate that uses only {count, S, Es0, Es1, CS}.  Also, at each solution, evaluate every certificate's
'kappa - (Psi-integral + P-term)' identity numerically: for a certificate using only CS, the Farkas inequality
int Psi dnu + pcol P <= kappa(S) < 0 must FAIL to apply (the solution's S lies outside the interval).
Usage: xn_posctl.py lin certfile nstarts seed
"""
import sys
import re
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import least_squares
from xn_feas import unpack, resid, margins

lin, cf, ns, seed = int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
k = 9
m = k - lin


def S_of_R(R):
    return (R + 0.5) ** 2 / 8 + (m - 1) / 32


certs = []
for line in open(cf):
    mm = re.match(r"\[(\d+), ([^,]+), ([^,]+), \[", line)
    Ra, Rb = float(Fr(mm.group(2))), float(Fr(mm.group(3)))
    kinds = set(re.findall(r'\["([A-Z0-9]+)"', line))
    certs.append((S_of_R(Ra), S_of_R(Rb), kinds))
rng = np.random.default_rng(seed)
found = 0
for it in range(ns):
    p0 = np.concatenate([rng.normal(0, 1.5, m), rng.normal(0, 1.0, m)])
    with np.errstate(all="ignore"):
        try:
            sol = least_squares(resid, p0, args=(k, lin, 2), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
        except Exception:
            continue
    if not np.all(np.isfinite(sol.fun)) or np.abs(sol.fun).max() > 1e-9:
        continue
    u, v = unpack(sol.x)
    mg, dm = margins(u, v)
    if mg <= 0 or dm < 1e-6:
        continue
    found += 1
    S = v.sum()
    hit = [c for c in certs if c[0] <= S <= c[1]]
    kinds = hit[0][2] if hit else None
    simple = kinds is not None and kinds <= {"CS"}
    print(f"N=2 solution {found}: S = {S:.6f}; certificate kinds there: {sorted(kinds) if kinds else 'NONE (outside cover)'}"
          f"  -> {'CONTRADICTION (bug!)' if simple else 'consistent'}", flush=True)
print(f"lin={lin}: {found} boundary solutions checked")
