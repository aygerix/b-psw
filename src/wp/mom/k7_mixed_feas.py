"""k = 7 families with m cubics (+ optional 2t+1) and invisible members of degree >= 7 (Theorem 11' hypotheses).
Feasibility of the four exact boundary equations (moment form, k7_mixed.gp) on the reality domain
8 v_i > (|u_i - 1/2| + 1/2)^2 (equivalent to c^2 < 8a and (a-c)^2 < 8a). Maximises the margin
delta = min_i [8 v_i - (|u_i - 1/2| + 1/2)^2] subject to the equations. delta_max < 0 => infeasible (numerically).
Usage: python3 k7_mixed_feas.py m lin starts"""
import sys, numpy as np
from scipy.optimize import minimize
m, lin, starts = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
rng = np.random.default_rng(7 + 10 * m + lin)
def mom(x):
    u, v = x[:m], x[m:]
    S = v.sum(); P1 = (v*u).sum(); P2 = (v*u**2).sum(); P3 = (v*u**3).sum(); M = (v**2).sum(); N = (v**2*u).sum()
    return S, P1, P2, P3, M, N
def eqs(x):
    S, P1, P2, P3, M, N = mom(x)
    if lin == 0:
        return np.array([0.5*S**2 - 4*S - 3*P1 - 2*P2 - P3 + 5.5*M + 4*N + 63/64,
                         6*S + 3*P1 + P2 - 2*M - 189/32,
                         0.5*S**2 - 10*S + 10*P1 - 5*P2 + P3 + 9.5*M - 4*N + 63/64,
                         10*S - 5*P1 + P2 - 2*M - 189/32])
    return np.array([0.5*S**2 - 2.5*S - 2*P1 - 1.5*P2 - P3 + 4.5*M + 4*N + 31/64,
                     4.5*S + 2.5*P1 + P2 - 2*M - 109/32,
                     0.5*S**2 - 7*S + 8*P1 - 4.5*P2 + P3 + 8.5*M - 4*N + 31/64,
                     8*S - 4.5*P1 + P2 - 2*M - 109/32])
# m = 1 special: pair terms vanish automatically since (S^2 - M)/2 = 0 for one member (moment form valid for all m >= 1)
def margins(x):
    u, v = x[:m], x[m:]
    return 8*v - (np.abs(u - 0.5) + 0.5)**2
best = -np.inf; bestx = None
for s in range(starts):
    v0 = rng.uniform(1/32, 1.0, m); u0 = 0.5 + rng.uniform(-1, 1, m) * np.sqrt(np.maximum(8*v0, 0.26)) * 0.5
    z0 = np.concatenate([u0, v0, [0.0]])
    cons = [{'type': 'eq', 'fun': lambda z: eqs(z[:-1])},
            {'type': 'ineq', 'fun': lambda z: margins(z[:-1]) - z[-1]}]
    try:
        r = minimize(lambda z: -z[-1], z0, constraints=cons, method='SLSQP', options={'maxiter': 2000, 'ftol': 1e-14})
    except Exception:
        continue
    z = r.x
    if np.abs(eqs(z[:-1])).max() < 1e-9:
        d = margins(z[:-1]).min()
        if d > best: best, bestx = d, z[:-1].copy()
print(f"m = {m}, lin = {lin}: best margin {best:.6e}" + ("" if bestx is None else
      f"  S = {bestx[m:].sum():.6f}  v = {np.round(np.sort(bestx[m:]), 6)}  u = {np.round(bestx[:m], 4)}"))
