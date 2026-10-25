"""Merge the FX lambda = 0.2 runs of the five unfinished E1 classes into run_bbc.sh format (for g80.py / summarize.py).
A split instance (sub-instances by fixed sides of rows 3, 5, 7, 11) gets M_ub = max over its sub-instances (each
configuration is feasible, with the same value, in the sub-instance matching its side pattern); it is reported only if
every sub-instance finished (DONE, rc 0).  Usage: python3 merge_l20.py L20DIR > l20_merged.out"""
import sys, os, glob, re
d = sys.argv[1]
POOL = {(-7, 2, 1): 4185, (-7, 1, 2): 4186, (-11, 2, 1): 3229, (-11, 2, 0): 3284, (-15, 2, 0): 6523}
SUBS = {(-15, 2, 0): ['FF', 'FL', 'LF_F', 'LF_L', 'LL_F_F', 'LL_F_L', 'LL_L'], (-11, 2, 0): ['F', 'L']}
def parse(fn):
    if not os.path.exists(fn): return None
    lam = [l for l in open(fn) if l.startswith('lam')]; rc = [l for l in open(fn) if l.startswith('rc')]
    if not lam or not rc or rc[-1].split()[1] != '0' or lam[-1].split()[2] != 'DONE': return None
    a = lam[-1].split(); return float(a[4]), int(a[12]), float(a[18]), lam[-1]
for key, n in POOL.items():
    D, j1, j2 = key
    base = os.path.join(d, f"D{D}_{j1}_{j2}")
    print(f"== {D} {j1} {j2} {n} 0 0 0")
    r = parse(base + '.out')
    if r: print(r[3].rstrip()); print("rc 0"); continue
    if key in SUBS:
        rs = [parse(f"{base}_{s}.out") for s in SUBS[key]]
        if all(rs):
            M = max(x[0] for x in rs); nodes = sum(x[1] for x in rs); t = sum(x[2] for x in rs)
            print(f"lam 0.200000 DONE M_ub {M:.6f} best {M:.6f} N 0 bits 0.000 nodes {nodes} exact 0 pruned 0 time {t:.2f} sb 0 "
                  f"sbclosures 0 inc 0 split {len(rs)}"); print("rc 0"); continue
    print("rc 1")
