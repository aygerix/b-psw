"""Independent recomputation of the certified small-b G_80 (adversarial FX check).
Parses every run file, keeps per (D, j1, j2) and lambda all DONE M_ub values (flags inconsistent duplicates), checks that
every instance of both index files has DONE at lambda 0.3, 0.25, 0.2, checks monotonicity M(l1) >= M(l2) for l1 < l2, and
recomputes max over b in [60, 275] of UB(b)/(b+1) with UB(b) = min(pool, min_lambda lambda b + M(lambda))."""
import sys, math, glob
from collections import defaultdict
runs = sys.argv[1:]
M = defaultdict(lambda: defaultdict(list)); hdrN = {}
for fn in runs:
    cur = None
    for l in open(fn):
        a = l.split()
        if not a: continue
        if a[0] == '==':
            cur = (int(a[1]), int(a[2]), int(a[3])); hdrN.setdefault(cur, set()).add(int(a[4]))
        elif a[0] == 'lam' and cur is not None and len(a) > 4 and a[2] == 'DONE':
            M[cur][round(float(a[1]), 6)].append((float(a[4]), fn.split('/')[-1]))
idx = {}
for fn, tag in (('in/inst_index.tsv', 'ten'), ('in/instxD_index.tsv', 'xD')):
    for l in open(fn):
        a = l.split('\t'); idx[(int(a[0]), int(a[1]), int(a[2]))] = (int(a[3]), tag)
print("instances in index files:", len(idx), " with run records:", len(set(M) & set(idx)))
missing = [(k, lam) for k in idx for lam in (0.3, 0.25, 0.2) if lam not in M.get(k, {})]
print("instance/lambda pairs WITHOUT a DONE record:", len(missing), missing[:10])
incons = [(k, lam, v) for k in M for lam, v in M[k].items() if max(x[0] for x in v) - min(x[0] for x in v) > 1e-6]
print("inconsistent duplicate DONE values:", len(incons), incons[:5])
hdrbad = [(k, hdrN[k], idx[k][0]) for k in idx if k in hdrN and hdrN[k] != {idx[k][0]}]
print("header pool size != index pool size:", hdrbad[:5])
mono = []
for k in M:
    ls = sorted(M[k]); vals = [min(x[0] for x in M[k][l]) for l in ls]
    for i in range(len(ls) - 1):
        if vals[i] + 1e-6 < vals[i + 1]: mono.append((k, ls[i], vals[i], ls[i + 1], vals[i + 1]))
print("monotonicity violations (M must not increase with lambda):", mono[:5], len(mono))
perD = defaultdict(lambda: (0, None)); allmax = (0, None)
for k, (n, tag) in idx.items():
    lv = {l: max(x[0] for x in v) for l, v in M.get(k, {}).items()}   # if duplicates, use the LARGEST (most conservative)
    if not lv: continue
    for b in range(60, 276):
        ub = min(n, min(l * b + m for l, m in lv.items()))
        r = ub / (b + 1)
        if r > perD[k[0]][0]: perD[k[0]] = (r, (k, n, b, round(ub, 4)))
        if r > allmax[0]: allmax = (r, (k, n, b, round(ub, 4)))
for D in sorted(perD, key=lambda d: (abs(d), d)):
    r, arg = perD[D]; print(f"  D {D:4d}: max r {r:.4f}  G80 >= {-math.log2(r):.3f}  binding {arg}")
print(f"OVERALL max r {allmax[0]:.4f}  certified G80 >= {-math.log2(allmax[0]):.4f}  at {allmax[1]}")
# without the lambda = 0.225 run
k = (-15, 2, 0); lv = {l: max(x[0] for x in v) for l, v in M[k].items() if abs(l - 0.225) > 1e-9}
rr = max((min(6523, min(l * b + m for l, m in lv.items())) / (b + 1), b) for b in range(60, 276))
print(f"D=-15 (2,0) without lambda=0.225: max r {rr[0]:.4f} at b={rr[1]} -> G80 {-math.log2(rr[0]):.3f}; M values {dict(sorted(lv.items()))}")
