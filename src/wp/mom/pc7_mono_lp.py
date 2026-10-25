"""Termwise certificate LP at fixed shapes: alpha with C0.alpha = 1, phi1_alpha >= 0 on (0, inf), Gam'_alpha >= 0 on (0, inf),
kap.alpha >= 0. (Sufficient for sum alpha_i B_i > 0 on all w > 0, any number of members.) Sampled positivity (log grid)."""
import numpy as np, sys
from scipy.optimize import linprog
from pc7_mono_lib import parts
X = np.exp(np.linspace(np.log(1e-5), np.log(1e5), 400))
def rows(r1, r2):
    C0, P1, A, kap = map(np.array, parts(r1, r2))
    R1 = np.array([[sum(P1[i][a] * x**(a - 1) for a in range(1, 10)) / (1 + x**8) for i in range(8)] for x in X])
    R2 = np.array([[sum((n + 1) * A[i][n] * x**n for n in range(5)) / (1 + x**4) for i in range(8)] for x in X])
    return C0, R1, R2, kap
def lp(r1, r2, supp=range(8)):
    C0, R1, R2, kap = rows(r1, r2); supp = list(supp)
    # vars alpha (len supp), delta; maximise delta s.t. R alpha >= delta, kap alpha >= 0, C0 alpha = 1, |alpha| <= 1e4
    M = np.vstack([R1[:, supp], R2[:, supp]]); n = len(supp)
    Aub = np.hstack([-M, np.ones((len(M), 1))]); bub = np.zeros(len(M))
    Aub = np.vstack([Aub, np.concatenate([-kap[supp], [0]])]); bub = np.append(bub, 0)
    Aeq = np.concatenate([C0[supp], [0]])[None, :]
    r = linprog(np.concatenate([np.zeros(n), [-1]]), A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=[1],
                bounds=[(-1e4, 1e4)] * n + [(None, 1)], method='highs')
    if r.status != 0: return None, None
    return r.x[-1], r.x[:n]
if __name__ == '__main__':
    for r1, r2 in [(-0.5, 0.5), (-0.3, 0.6), (0.1, 0.5), (-0.9, 0.2), (0.3, 0.9), (-0.8, -0.1), (0.0, 0.2), (0.6, 0.7)]:
        d, a = lp(r1, r2)
        print(f"r = ({r1:+.2f},{r2:+.2f}) margin {d}  alpha {None if a is None else np.round(a, 4)}", flush=True)
