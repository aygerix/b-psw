"""Univariate reduction: with c_n (n = 2..6) the Hankel coefficients of b_alpha (b = sum_n c_n sum_{a+b=n} w^a v^b),
Chebyshev/positivity bounds p_a p_b <= 7 p_{a+b} and p_a p_b >= p_{a+b} give
   E_alpha >= kap + sum_j phi(w_j),  phi(w) = a_alpha(w) - 3 sum_n (n-1) s_n w^n,  s_n >= max(0, -c_n),
so phi >= -kap/7 on (0, inf) (and lam3 >= 0) certifies the shape.  LP in (alpha, s).  Usage: python3 pc7_cert_uni.py ngrid"""
import sys, numpy as np
from scipy.optimize import linprog
from pc7_cert_lib import load, at_shape, LAM3
D = load()
WG = np.exp(np.linspace(np.log(1e-6), np.log(1e6), 361))
def test(R1, R2):
    kap, Ac, Bc = at_shape(D, R1, R2)
    C = np.array([[Bc[i, 1, n - 1] for n in range(2, 7)] for i in range(8)])     # (8, 5): c_{i,n}
    Ra = np.array([[np.polyval(Ac[i][::-1], x) for i in range(8)] for x in WG])
    Pn = np.array([[-3 * (n - 1) * x ** n for n in range(2, 7)] for x in WG])
    sc = np.maximum(np.maximum(np.abs(Ra).max(axis=1), np.abs(Pn).max(axis=1)), 1)[:, None]
    # vars: alpha (8), s (5).  phi(w) + kap/7 >= 0;  s_n + c_n(alpha) >= 0;  s >= 0;  lam3 >= 0
    A1 = -np.hstack([Ra + kap[None, :] / 7, Pn]) / sc
    A2 = -np.hstack([C.T, np.eye(5)])
    A3 = -np.concatenate([LAM3, np.zeros(5)])[None, :]
    Aub = np.vstack([A1, A2, A3])
    r = linprog(np.zeros(13), A_ub=Aub, b_ub=np.zeros(len(Aub)), A_eq=np.concatenate([kap, np.zeros(5)])[None, :], b_eq=[1.0],
                bounds=[(-1e5, 1e5)] * 8 + [(0, 1e7)] * 5, method='highs')
    return r.status == 0
ng = int(sys.argv[1]); grid = np.linspace(-0.9, 0.9, ng)
shapes = [(a, b) for a in grid for b in grid if b >= a - 1e-12]
res = [(s, test(*s)) for s in shapes]
bad = [(round(s[0], 2), round(s[1], 2)) for s, ok in res if not ok]
print(f"uni: {len(res) - len(bad)}/{len(res)} shapes certified; failing {bad}")
