"""Constant alpha for the reality-restricted problem at one shape (max margin over 7-member configurations in [lo, hi]^7),
written as a rational vector for the exact branch-and-bound pc7_cert_bb.gp.  Usage: python3 pc7_cert_ralpha.py r1 r2 lo hi nconf out"""
import sys, numpy as np
from fractions import Fraction as Fr
from scipy.optimize import linprog
from pc7_cert_dlib import load, Evec
at = load(); rng = np.random.default_rng(17)
r1, r2 = float(Fr(sys.argv[1])), float(Fr(sys.argv[2])); lo, hi = float(Fr(sys.argv[3])), float(Fr(sys.argv[4])); nconf = int(sys.argv[5])
K, A, C = at(r1, r2)
CF = [rng.uniform(lo, hi, 7) for _ in range(nconf)]
for _ in range(nconf):
    m = rng.integers(1, 4); vals = rng.uniform(lo, hi, m); idx = rng.integers(0, m, 7); CF.append(vals[idx])
grid = np.linspace(lo, hi, 9)
for x in grid:                                   # corners / edges of the box
    for y in grid: CF.append(np.array([x] * 3 + [y] * 4)); CF.append(np.array([x] * 6 + [y]))
D = np.array([Evec(K, A, C, c) for c in CF]); Dn = D / np.abs(D).sum(axis=1, keepdims=True)
for it in range(3):
    r = linprog(np.concatenate([np.zeros(8), [-1]]), A_ub=np.hstack([-Dn, np.ones((len(Dn), 1))]), b_ub=np.zeros(len(Dn)), bounds=[(-1, 1)] * 8 + [(None, None)], method='highs')
    al = r.x[:8]; t = r.x[-1]
    # add local minimisers of the normalised E_alpha
    from scipy.optimize import minimize
    f = lambda z: (lambda E: E @ al / np.abs(E).sum())(Evec(K, A, C, lo + (hi - lo) / (1 + np.exp(-z))))
    worst = np.argsort(Dn @ al)[:10]
    for j in worst:
        z0 = -np.log((hi - lo) / (CF[j] - lo + 1e-12) - 1 + 1e-12)
        res = minimize(lambda z: f(np.clip(z, -30, 30)), np.clip(np.nan_to_num(z0), -20, 20), method='Nelder-Mead', options={'maxiter': 3000})
        w = lo + (hi - lo) / (1 + np.exp(-np.clip(res.x, -30, 30))); Ew = Evec(K, A, C, w)
        if np.all(np.isfinite(Ew)) and np.abs(Ew).sum() > 0: CF.append(w); D = np.vstack([D, Ew])
    Dn = D / np.abs(D).sum(axis=1, keepdims=True)
    print(f"iter {it}: margin {t:+.4e}, min over local searches {min(f(-np.log((hi - lo) / (c - lo + 1e-12) - 1 + 1e-12)) for c in CF[-10:]):+.4e}", flush=True)
alq = [Fr(x).limit_denominator(10**6) for x in al]
with open(sys.argv[6], 'w') as fo: fo.write("ALPHA = [" + ", ".join(str(x) for x in alq) + "];\n")
print("alpha =", [str(x) for x in alq])
