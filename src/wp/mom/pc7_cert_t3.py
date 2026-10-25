import sys, itertools, numpy as np
from fractions import Fraction as Fr
from scipy.optimize import linprog
from pc7_cert_exact import *
rng = np.random.default_rng(3)
G = grid(-14, 14)
for S in load('pc7_cert_t3.dat'):
    K = np.array([float(x) for x in S['K']])
    R = nrm(pair_rows(S, G)); rp = lp_feas(R, K)
    # full: random 7-member configurations, exact
    rows = []
    for _ in range(600):
        ws = [Fr(2)**int(e) * Fr(int(m), 64) for e, m in zip(rng.integers(-14, 14, 7), rng.integers(64, 128, 7))]
        E = [S['K'][i] + sum(aval(S, i, x) for x in ws) + sum(bval(S, i, x, y) for x, y in itertools.combinations(ws, 2))
             + LAM3[i] * sum(x * y * z for x, y, z in itertools.combinations(ws, 3)) for i in range(8)]
        rows.append(E)
    Rf = nrm(rows)
    rf = linprog(np.zeros(8), A_ub=-Rf, b_ub=-1e-9 * np.ones(len(Rf)), bounds=[(-1e6, 1e6)] * 8, method='highs')
    print(f"r = ({float(S['r'][0]):+.2f},{float(S['r'][1]):+.2f}): pair {'OK' if rp.status == 0 else 'FAIL'}   full {'OK' if rf.status == 0 else 'FAIL'}", flush=True)
