#!/usr/bin/env python3
"""src/cs/e3/check_scan.py -- independent check of an E3 scan file with the exact verifier src/bpsw.py.
Each line: n p1 p2 p3 rho d X D s2 sl kor dF dL1 dL2 dL3 dLt  (from rho_scan.gp / p23_scan.gp).
Checks n = p1 p2 p3, recomputes strong_fermat2, strong_lucas, Selfridge D with bpsw.py and compares with gp.
Reports the Lucas-defect statistics and the heuristic expected hit count sum 1/dLt over all-inert tuples.
Usage: check_scan.py SCANFILE [path/to/src]"""
import sys, math
sys.path.insert(0, sys.argv[2] if len(sys.argv) > 2 else '../..')
import bpsw
rows = [l.split() for l in open(sys.argv[1]) if l.strip()]
mism = hits = 0; inert = 0; exp_hits = 0.0; dl1_one_inert = 0; lr = []
for r in rows:
    n, p1, p2, p3 = map(int, r[:4]); D = int(r[7]); s2, sl = int(r[8]), int(r[9])
    dL1, dL2, dL3, dLt = map(int, r[12:16])
    assert n == p1 * p2 * p3
    F2 = bpsw.strong_fermat2(n); SL = bpsw.strong_lucas(n); Dp = bpsw.selfridge_D(n)
    if (F2, SL, Dp) != (bool(s2), bool(sl), D):
        mism += 1
    if F2 and SL:
        hits += 1; print("BPSW PASS (bpsw.py):", n, p1, p2, p3)
    if all(bpsw.jacobi(D % p, p) == -1 for p in (p1, p2, p3)):
        inert += 1; exp_hits += 1.0 / dLt
        dl1_one_inert += (dL1 == 1)
        lr.append(math.log(dLt) / math.log(n) if dLt > 1 else 0.0)
lr.sort()
print(f"rows {len(rows)}  gp/bpsw.py mismatches {mism}  bpsw.py passes {hits}")
print(f"all-inert tuples {inert}; of these member-1 Lucas defect 1: {dl1_one_inert}")
if lr:
    print(f"inert: log dL/log n  min {lr[0]:.3f}  median {lr[len(lr)//2]:.3f}  max {lr[-1]:.3f}")
    print(f"heuristic expected strong-Lucas absorptions sum_inert 1/dL = {exp_hits:.3e}")
