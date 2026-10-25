#!/usr/bin/env python3
"""ns line: summary table over ns_out/v3_*.npz.  For each run: starts, planted recall, and for NON-DEGENERATE endpoints
(member separation, special-point gaps (incl. conjugate pairs = distance to the reality boundary), min|a| all >= DL)
the min / 1% / 10% / median of the identity residual, the remainder residual and their max.
Usage: ns_summary.py [DL=0.05] [glob]
"""
import sys, glob
import numpy as np

DL = float(sys.argv[1]) if len(sys.argv) > 1 else 0.05
pat = sys.argv[2] if len(sys.argv) > 2 else 'ns_out/v3_*.npz'
AL = float(sys.argv[3]) if len(sys.argv) > 3 else DL     # separate threshold for min|a| (near-linear collapse)
print(f"non-degenerate := sep >= {DL}, gap >= {DL}, min|a| >= {AL}")
print(f"{'run':45s} {'n':>6s} {'nd':>6s} {'recall':>8s} | {'min both':>9s} {'1%':>9s} {'10%':>9s} {'50%':>9s} | "
      f"{'min id':>9s} {'min rem':>9s} | {'all: min both':>13s}")
for f in sorted(glob.glob(pat)):
    d = np.load(f)
    rr, rem, sep, gap, mina = d['rr'], d['rem'], d['sep'], d['gap'], d['mina']
    both = np.maximum(rr, np.nan_to_num(rem, nan=np.inf))
    nd = (sep >= DL) & (gap >= DL) & (mina >= AL)
    planted = '_PL' in f
    rec = ''
    if planted:
        rec = f"{np.mean((rr < 1e-10) & (d['dstar'] < 1e-6)):.3f}"
    name = f.split('/')[-1].replace('.npz', '')
    if nd.any():
        b = both[nd]
        q = np.quantile(b, [0, 0.01, 0.1, 0.5])
        print(f"{name:45s} {len(rr):6d} {nd.sum():6d} {rec:>8s} | " + " ".join(f"{z:9.2e}" for z in q) +
              f" | {rr[nd].min():9.2e} {np.nanmin(rem[nd]) if np.isfinite(rem[nd]).any() else np.nan:9.2e} | {both.min():13.2e}")
    else:
        print(f"{name:45s} {len(rr):6d} {0:6d} {rec:>8s} | (no non-degenerate endpoint) | all: min both {both.min():.2e}")
