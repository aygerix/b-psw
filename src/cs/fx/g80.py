"""Certified small-b G_80 per D from bbc / bbc_fx outputs (run_bbc.sh / run_fx.sh format; later files override earlier
ones per (instance, lambda)).  UB_c(b) = min(pool size, min over finished lambda of lambda*b + M_ub(lambda)) (weak duality, every b);
ratio r_c(b) = UB_c(b)/(b+1); G_80(D) >= -log2 max_{c, 60 <= b <= 275} r_c(b)  (existence regime N_req = b + 1 for b <= 275).
Also prints the floor version (N_max is an integer, so N_max(b) <= floor(UB(b))), the binding instance, and instances with
no finished lambda (not certified).  Usage: python3 g80.py OUT1 [OUT2 ...] [--extra KEY=OUTFILE ...]"""
import sys, math
from collections import defaultdict
inst = {}
for fn in sys.argv[1:]:
    cur = None
    for l in open(fn):
        a = l.split()
        if not a: continue
        if a[0] == '==':
            cur = (int(a[1]), int(a[2]), int(a[3])); inst.setdefault(cur, dict(n=int(a[4]), lam={}))
        elif a[0] == 'lam' and cur and a[2] == 'DONE':
            inst[cur]['lam'][float(a[1])] = float(a[4])
B = range(60, 276)
perD = defaultdict(lambda: dict(r=0, rf=0, arg=None, argf=None, nocert=[], ninst=0, lams=defaultdict(int)))
for key, v in inst.items():
    d = perD[key[0]]; d['ninst'] += 1
    if not v['lam']: d['nocert'].append((key[1:], v['n'])); continue
    for lam in v['lam']: d['lams'][lam] += 1
    for b in B:
        ub = min(min(lam * b + M for lam, M in v['lam'].items()), v['n'])   # trivial bound N <= pool size
        r = ub / (b + 1); rf = math.floor(ub + 1e-9) / (b + 1)
        if r > d['r']: d['r'] = r; d['arg'] = (key[1:], v['n'], b, round(ub, 3), min(v['lam'], key=lambda L: L * b + v['lam'][L]) if ub < v['n'] else 'pool')
        if rf > d['rf']: d['rf'] = rf; d['argf'] = (key[1:], b, math.floor(ub + 1e-9))
print(f"{'D':>5} {'inst':>4}  {'max r':>7} {'G80>=':>6}  {'floor r':>7} {'G80f>=':>6}  binding (class, pool, b, UB, lambda)   lambdas finished   not certified")
allr = 0
for D in sorted(perD, key=lambda x: (abs(x), x)):
    d = perD[D]; allr = max(allr, d['r'])
    g = -math.log2(d['r']) if d['r'] > 0 else float('inf'); gf = -math.log2(d['rf']) if d['rf'] > 0 else float('inf')
    print(f"{D:5d} {d['ninst']:4d}  {d['r']:.4f} {g:6.3f}  {d['rf']:.4f} {gf:6.3f}  {d['arg']}   "
          f"{dict(sorted(d['lams'].items(), reverse=True))}   {d['nocert'][:6]}")
print(f"overall: max r = {allr:.4f}, certified G80 >= {-math.log2(allr):.3f}")
