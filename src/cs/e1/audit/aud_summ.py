"""Independent parse of the bbc pass outputs: per instance M_ub(lambda) for finished lambdas; maxima over admissible /
all classes; unfinished instances; certified envelope UB(b) = min_lambda (lambda b + M) and worst UB(b)/(b+1);
D = 5 classes.  Usage: python3 aud_summ.py passfiles..."""
import sys
from collections import defaultdict
inst = {}; cur = None
for fn in sys.argv[1:]:
    for l in open(fn):
        a = l.split()
        if not a: continue
        if a[0] == '==': cur = (int(a[1]), int(a[2]), int(a[3])); inst.setdefault(cur, {'n': int(a[4]), 'lam': {}, 'rc': []})
        elif a[0] == 'lam': inst[cur]['lam'][float(a[1])] = (a[2], float(a[4]) if a[2] == 'DONE' else None, int(a[8]), float(a[10]))
        elif a[0] == 'rc': inst[cur]['rc'].append(int(a[1]))
adm = lambda k: min(k[1], k[2]) <= 1
print(f"instances {len(inst)}, admissible {sum(1 for k in inst if adm(k))}")
for lam in (0.3, 0.25, 0.2):
    done = {k: v['lam'][lam][1] for k, v in inst.items() if lam in v['lam'] and v['lam'][lam][0] == 'DONE'}
    ma = max((M, k) for k, M in done.items() if adm(k)); mall = max((M, k) for k, M in done.items())
    nd = [k for k in inst if k not in done]
    print(f"lambda {lam}: done {len(done)} (admissible {sum(1 for k in done if adm(k))}); max M admissible {ma}, all {mall}; "
          f"not done: {[(k, inst[k]['n']) for k in nd]}")
for b in (60, 100, 140, 214, 227, 275):
    w = max(((min(lam * b + v['lam'][lam][1] for lam in v['lam'] if v['lam'][lam][0] == 'DONE')) / (b + 1), k)
            for k, v in inst.items() if adm(k))
    wall = max(((min(lam * b + v['lam'][lam][1] for lam in v['lam'] if v['lam'][lam][0] == 'DONE')) / (b + 1), k)
            for k, v in inst.items())
    print(f"b = {b}: worst UB/(b+1) admissible {w[0]:.4f} at {w[1]}; all {wall[0]:.4f} at {wall[1]}")
print("D = 5 classes, M(0.2):", sorted(((v['lam'].get(0.2, (None, None))[1], k) for k, v in inst.items() if k[0] == 5),
                                   key=lambda t: -(t[0] or 0))[:6])
print("Q = +-2^k D (-7, -15, 17, -31) vs others: max M(0.25) per D:",
      {D: max(v['lam'][0.25][1] for k, v in inst.items() if k[0] == D and adm(k)) for D in (5, -7, -11, 13, -15, 17, -19, -23, 29, -31)})
print("max M(0.2) per D (admissible, done):",
      {D: max((v['lam'][0.2][1] for k, v in inst.items() if k[0] == D and adm(k) and 0.2 in v['lam']), default=None) for D in (5, -7, -11, 13, -15, 17, -19, -23, 29, -31)})
