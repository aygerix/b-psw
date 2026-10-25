"""Termwise LP restricted to the REALITY interval I = I- cap I+ of w (K-_j, K+_j without real roots):
 q+(w) = lam^2 w^2 - 4 (mu + lam e1) w + 4 < 0,  q-(w) = lam^2 w^2 - (2 lam e1 + 4 mu) w + (r1 - r2)^2 < 0."""
import numpy as np, sys
from scipy.optimize import linprog
from pc7_mono_lib import parts, shape_env
def interval(r1, r2):
    d = shape_env(r1, r2); lam, mu, e1 = d['lam'], d['mu'], d['e1']
    def roots(a, b, c):
        if abs(a) < 1e-14: return (c / -b if b < 0 else -np.inf, np.inf) if True else None
        D = b * b - 4 * a * c
        if D <= 0: return None
        s = np.sqrt(D); return ((-b - s) / (2 * a), (-b + s) / (2 * a))
    Ip = roots(lam**2, -4 * (mu + lam * e1), 4.0); Im = roots(lam**2, -(2 * lam * e1 + 4 * mu), (r1 - r2)**2)
    if Ip is None or Im is None: return None
    lo, hi = max(Ip[0], Im[0]), min(Ip[1], Im[1])
    return (lo, hi) if lo < hi else None
def lp(r1, r2, supp=range(8), ns=300, verbose=False):
    I = interval(r1, r2)
    if I is None: return 'EMPTY', None, None
    lo, hi = I; hi_s = min(hi, lo * 1e6)
    X = np.exp(np.linspace(np.log(lo), np.log(hi_s), ns))
    C0, P1, A, kap = map(np.array, parts(r1, r2)); supp = list(supp); n = len(supp)
    R1 = np.array([[sum(P1[i][a] * x**(a - 1) for a in range(1, 10)) / (1 + x**8) for i in supp] for x in X])
    R2 = np.array([[sum((m + 1) * A[i][m] * x**m for m in range(5)) / (1 + x**4) for i in supp] for x in X])
    M = np.vstack([R1, R2])
    Aub = np.hstack([-M, np.ones((len(M), 1))]); bub = np.zeros(len(M))
    Aub = np.vstack([Aub, np.concatenate([-kap[supp], [0]])]); bub = np.append(bub, 0)
    Aeq = np.concatenate([C0[supp], [0]])[None, :]
    r = linprog(np.concatenate([np.zeros(n), [-1]]), A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=[1],
                bounds=[(-1e4, 1e4)] * n + [(None, 1)], method='highs')
    if r.status != 0: return 'LPFAIL', I, None
    return r.x[-1], I, r.x[:n]
if __name__ == '__main__':
    np.set_printoptions(linewidth=200, precision=4)
    for r1, r2 in [(-0.3, 0.6), (0.1, 0.5), (0.0, 0.2), (-0.9, 0.2), (0.3, 0.9), (-0.8, -0.1), (0.6, 0.7), (-0.5, 0.5), (-0.1, 0.1), (0.85, 0.95), (-0.95, 0.9)]:
        d, I, a = lp(r1, r2)
        print(f"r = ({r1:+.2f},{r2:+.2f}) I = {None if I is None else np.round(I, 5)} margin {d}  alpha {a}", flush=True)
        for s in [[0], [3], [0, 3], [0, 1, 2, 3]]:
            d2, _, a2 = lp(r1, r2, supp=s)
            print(f"      supp {[i + 1 for i in s]}: margin {d2}  alpha {a2}")
