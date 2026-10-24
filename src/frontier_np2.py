"""Refined frontier: finer lambda grid, more iterations, more seeds (reuses frontier_np).
Usage: python3 frontier_np2.py z pairs_file iters seeds class_rank"""
import sys, math, json
from collections import defaultdict
import numpy as np
from sympy import primerange
import frontier_np as F
from pool import build
z = int(sys.argv[1]); rows = build(z, sys.argv[2]); iters = int(sys.argv[3]); seeds = int(sys.argv[4]); rk = int(sys.argv[5])
qs = [q for q in primerange(3, z + 1) if q != 5]
lq1 = np.array([math.log2(q - 1) for q in qs]); lq = np.array([math.log2(q) for q in qs])
bycls = defaultdict(list)
for r in rows: bycls[r['cls']].append(r)
c = sorted(bycls, key=lambda c: -len(bycls[c]))[rk]; R = bycls[c]; P = F.prep(R, qs)
F.qs_p = {t: R[t]['p'] for t in P[4]}
out = []
for lam in (0.0, 0.03, 0.06, 0.08, 0.1, 0.12, 0.14, 0.17):
    b = max((F.anneal_l(P, qs, lq1, lq, lam, iters, s) for s in range(seeds)), key=lambda t: t[0])
    A = [q for q in qs if b[1][q] == 'A']; B = [q for q in qs if b[1][q] == 'B']
    print(f"z={z} class {c} lam {lam}: log2|G| {b[3]:.1f} N {b[2]} ratio {b[2]/b[3]:.3f}  A={A} B={B}", flush=True)
    out.append((lam, round(b[3], 1), b[2], A, B))
json.dump(out, open(f"../data/frontier2_z{z}_c{rk}.json", "w"))
