"""Reality-restricted pair / triple relaxations with EXACT rows (Fractions) on a rational superset [lo, hi] of the reality
interval I(shape).  Prints max margins and the pair-optimal alpha (as floats, for rationalisation).
Usage: python3 pc7_cert_rpair.py dumpfile r1 r2 lo hi n2 n3"""
import sys, itertools, numpy as np
from fractions import Fraction as Fr
from scipy.optimize import linprog
from pc7_cert_exact import load, aval, bval, LAM3
S = [s for s in load(sys.argv[1]) if s['r'] == (Fr(sys.argv[2]), Fr(sys.argv[3]))][0]
lo, hi = Fr(sys.argv[4]), Fr(sys.argv[5]); n2, n3 = int(sys.argv[6]), int(sys.argv[7])
G2 = [lo + (hi - lo) * Fr(k, n2) for k in range(n2 + 1)]; G3 = [lo + (hi - lo) * Fr(k, n3) for k in range(n3 + 1)]
kap = np.array([float(x) for x in S['K']])
def mm(rows, need_lam):
    R = np.array([[float(x) for x in r] for r in rows]); Rn = R / np.abs(R).sum(axis=1, keepdims=True); L3 = np.array(LAM3, float)
    Aub = [np.hstack([-Rn, Rn, np.ones((len(Rn), 1))]), np.concatenate([-kap, kap, [0]])[None, :], np.concatenate([np.ones(16), [0]])[None, :]]
    if need_lam: Aub.insert(1, np.concatenate([-L3, L3, [0]])[None, :])
    Aub = np.vstack(Aub); bub = np.zeros(len(Aub)); bub[-1] = 1
    r = linprog(np.concatenate([np.zeros(16), [-1]]), A_ub=Aub, b_ub=bub, bounds=[(0, None)] * 16 + [(None, None)], method='highs')
    return r.x[-1], r.x[:8] - r.x[8:16]
av = {x: [aval(S, i, x) for i in range(8)] for x in set(G2) | set(G3)}
prow = [[S['K'][i] / 21 + (av[x][i] + av[y][i]) / 6 + bval(S, i, x, y) for i in range(8)] for x, y in itertools.combinations_with_replacement(G2, 2)]
tp, ap = mm(prow, True)
print(f"pair margin {tp:+.4e}  alpha = {list(np.round(ap / np.abs(ap).max(), 8))}", flush=True)
bcache = {}
def bb(x, y):
    if (x, y) not in bcache: bcache[(x, y)] = [bval(S, i, x, y) for i in range(8)]
    return bcache[(x, y)]
trow = [[S['K'][i] / 35 + (av[x][i] + av[y][i] + av[z][i]) / 15 + (bb(x, y)[i] + bb(x, z)[i] + bb(y, z)[i]) / 5 + LAM3[i] * x * y * z for i in range(8)]
        for x, y, z in itertools.combinations_with_replacement(G3, 3)]
tt, at_ = mm(trow, False)
print(f"triple margin {tt:+.4e}  alpha = {list(np.round(at_ / np.abs(at_).max(), 8))}")
