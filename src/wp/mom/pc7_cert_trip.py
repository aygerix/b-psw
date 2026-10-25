"""Exact-row LPs (D-basis) at the shapes of a dump file: pair and triple relaxations, with max margin.
  pair  : pi(w, v) = kap/21 + (a(w) + a(v))/6 + b(w, v)
  trip  : tau(u, v, w) = kap/35 + (a(u) + a(v) + a(w))/15 + (b(u, v) + b(u, w) + b(v, w))/5 + lam3 u v w
Both sufficient for E_alpha > 0 on (0, inf)^7 when >= 0 everywhere with kap_alpha > 0 (pair also needs lam3 >= 0).
Usage: python3 pc7_cert_trip.py dumpfile glo ghi [step3]"""
import sys, itertools, numpy as np
from fractions import Fraction as Fr
from scipy.optimize import linprog
from pc7_cert_exact import load, aval, bval, LAM3
def gridf(lo, hi, step=1):
    g = []
    for e in range(lo, hi + 1, step): g += [Fr(2)**e, Fr(3, 2) * Fr(2)**e]
    return sorted(set(g))
def maxmargin(rows, kap, need_lam):
    R = np.array([[float(x) for x in r] for r in rows]); N = np.abs(R).sum(axis=1, keepdims=True); Rn = R / N
    L3 = np.array(LAM3, float)
    Aub = [np.hstack([-Rn, Rn, np.ones((len(Rn), 1))])]
    if need_lam: Aub.append(np.concatenate([-L3, L3, [0]])[None, :])
    Aub.append(np.concatenate([-kap, kap, [0]])[None, :]); Aub.append(np.concatenate([np.ones(16), [0]])[None, :])
    Aub = np.vstack(Aub); bub = np.zeros(len(Aub)); bub[-1] = 1
    r = linprog(np.concatenate([np.zeros(16), [-1]]), A_ub=Aub, b_ub=bub, bounds=[(0, None)] * 16 + [(None, None)], method='highs')
    if r.status != 0: return None, None
    al = r.x[:8] - r.x[8:16]; return r.x[16], al / np.abs(al).max()
lo, hi = int(sys.argv[2]), int(sys.argv[3]); step3 = int(sys.argv[4]) if len(sys.argv) > 4 else 2
G2 = gridf(lo, hi); G3 = gridf(lo, hi, step3)
for S in load(sys.argv[1]):
    kap = np.array([float(x) for x in S['K']])
    av = {x: [aval(S, i, x) for i in range(8)] for x in G2}
    bv = {}
    def bb(x, y):
        if (x, y) not in bv: bv[(x, y)] = [bval(S, i, x, y) for i in range(8)]
        return bv[(x, y)]
    prow = [[S['K'][i] / 21 + (av[x][i] + av[y][i]) / 6 + bb(x, y)[i] for i in range(8)] for x, y in itertools.combinations_with_replacement(G2, 2)]
    tp, ap = maxmargin(prow, kap, True)
    trow = [[S['K'][i] / 35 + (av[x][i] + av[y][i] + av[z][i]) / 15 + (bb(x, y)[i] + bb(x, z)[i] + bb(y, z)[i]) / 5 + LAM3[i] * x * y * z for i in range(8)]
            for x, y, z in itertools.combinations_with_replacement(G3, 3)]
    tt, at_ = maxmargin(trow, kap, False)
    f = lambda t: 'infeasible' if t is None else f"{t:+.3e}"
    print(f"r = ({float(S['r'][0]):+.3f},{float(S['r'][1]):+.3f}): pair margin {f(tp)}   triple margin {f(tt)}  ({len(prow)}/{len(trow)} rows)", flush=True)
