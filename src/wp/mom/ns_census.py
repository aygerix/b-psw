#!/usr/bin/env python3
"""ns line: census of REAL "-" half solutions inside the reality domain D (g_i AND h_i without real roots), from one or
more half runs (ns_run3.py k lin real - ...).  A real Williams family is exactly a half solution x in D whose
reflection sigma(x) = (a, a - c) is again a half solution (the "+" half of x is the "-" half of sigma(x)).
For each distinct half solution: multiplicity (how many starts found it), its "+" residual (identity and per-member
remainder), the reality margin, and whether sigma(x) is in the found list.  Chao1 estimate of the total number.
Usage: ns_census.py k lin file1.npz [file2.npz ...]
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ns_core as nc
import ns_core2 as n2
import ns_core3 as n3

k, lin = int(sys.argv[1]), int(sys.argv[2])
files = [f for f in sys.argv[3:] if not f.startswith('--')]
NOMARGIN = '--nomargin' in sys.argv          # count real half solutions anywhere (not only inside D)
A = []; C = []; RR = []
for f in files:
    d = np.load(f)
    A.append(d['a']); C.append(d['c']); RR.append(d['rr'])
    print(f"{f}: {len(d['rr'])} runs")
a = np.concatenate(A); c = np.concatenate(C); rr = np.concatenate(RR)
m = a.shape[1]
y = n3.eval_points(k)
with np.errstate(all='ignore'):
    rem_m = np.abs(n2.rem_resid(a, c, lin, '-').real).max(axis=1)
    gap = n2.all_gap(a, c, '-+')[1]
sep = nc.pair_sep(a, c, lin)
mg = nc.margins(a, c)
ok = (rr < 1e-9) & (rem_m < 1e-7) & (sep > 1e-4) & ((mg > 1e-6) | NOMARGIN) & (np.abs(a).min(axis=1) > 1e-3)
if NOMARGIN:
    print("NOTE: --nomargin: solutions outside the reality domain D are counted too")
print(f"total runs {len(rr)}; converged half solutions in D (identity < 1e-9, remainder < 1e-7, sep > 1e-4, margin > 1e-6): "
      f"{ok.sum()}")
print(f"   (identity < 1e-9 but rejected: {((rr < 1e-9) & ~ok).sum()}; of those margin <= 1e-6: {((rr < 1e-9) & ~ok & (mg <= 1e-6)).sum()},"
      f" sep <= 1e-4: {((rr < 1e-9) & ~ok & (sep <= 1e-4)).sum()})")
idx = np.nonzero(ok)[0]


def key(ai, ci):
    pts = sorted(zip(np.round(ai.real, 5), np.round(ci.real, 5)))
    return tuple(pts)


groups = {}
for i in idx:
    kk = key(a[i], c[i])
    # tolerant matching: compare with existing keys
    found = None
    for g in groups:
        if max(abs(p[0] - q[0]) + abs(p[1] - q[1]) for p, q in zip(g, kk)) < 1e-4 * (1 + max(abs(p[0]) for p in kk)):
            found = g; break
    if found is None:
        groups[kk] = [i]
    else:
        groups[found].append(i)
mult = np.array([len(v) for v in groups.values()])
print(f"distinct half solutions (up to permutation): {len(groups)}; multiplicity histogram {{mult: count}}: "
      f"{dict(zip(*np.unique(mult, return_counts=True))) if len(mult) else {}}")
if len(mult):
    f1 = (mult == 1).sum(); f2 = (mult == 2).sum()
    chao = len(groups) + (f1 * f1 / (2 * f2) if f2 > 0 else f1 * (f1 - 1) / 2)
    print(f"Chao1 estimate of the number of distinct half solutions reachable by this sampler: {chao:.1f}")
keys = list(groups.keys())
reps = np.array([groups[g][0] for g in keys], dtype=int)
if len(reps):
    ar, cr = a[reps], c[reps]
    with np.errstate(all='ignore'):
        idp = np.abs(n3.resid(ar, cr, lin, y, '+', True)).max(axis=1)
        Rp = np.abs(n2.rem_resid(ar, cr, lin, '+').real)                       # (n, 2m): r0 (m), r1 (m)
        remp = Rp.max(axis=1)
        mm = ar.shape[1]
        remp_member = np.maximum(Rp[:, :mm], Rp[:, mm:])                         # per-member '+' residual
        remp_min = remp_member.min(axis=1)
        ids = np.abs(n3.resid(ar, ar - cr, lin, y, '-', True)).max(axis=1)   # '-' residual of the reflection
    order = np.argsort(remp)
    print(f"'+' residual over the distinct half solutions: min identity {idp.min():.3e}, min remainder {remp.min():.3e}")
    print(f"   quantiles of '+' remainder residual (0, 10, 50, 90%): {np.quantile(remp, [0, .1, .5, .9])}")
    print(f"   best single member '+' residual (partial family, J = 1): min {remp_min.min():.3e}, quantiles (0,10,50%) "
          f"{np.quantile(remp_min, [0, .1, .5])}")
    for j in order[:min(10, len(order))]:
        sk = key(ar[j], ar[j] - cr[j])
        refl_found = any(max(abs(p[0] - q[0]) + abs(p[1] - q[1]) for p, q in zip(g, sk)) < 1e-4 * (1 + max(abs(p[0]) for p in sk)) for g in keys)
        print(f"  #{j}: mult {len(groups[keys[j]])}, '+' identity {idp[j]:.3e}, '+' remainder {remp[j]:.3e}, "
              f"reflection '-' identity {ids[j]:.3e}, reflection in list: {refl_found}, margin {mg[reps[j]]:.3e}, sep {sep[reps[j]]:.3e}")
        print(f"      a = {np.array2string(ar[j].real, precision=9, max_line_width=400)}")
        print(f"      c = {np.array2string(cr[j].real, precision=9, max_line_width=400)}")
    np.savez_compressed(f"ns_out/census_k{k}_l{lin}" + ("_all" if NOMARGIN else "") + ".npz", a=ar, c=cr, mult=np.array([len(groups[g]) for g in keys]),
                        idp=idp, remp=remp, remp_member=remp_member)
