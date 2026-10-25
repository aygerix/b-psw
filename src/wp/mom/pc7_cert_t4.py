"""At a shape where the pair relaxation fails: (1) which (w, v) pairs obstruct it (dual of the pair LP);
(2) full-LP alpha (max margin over random + structured 7-member configurations) and the minimising configurations of E_alpha."""
import sys, itertools, numpy as np
from fractions import Fraction as Fr
from scipy.optimize import linprog, minimize
from pc7_cert_exact import *
S = [s for s in load(sys.argv[1]) if s['r'] == (Fr(sys.argv[2]), Fr(sys.argv[3]))][0]
K = np.array([float(x) for x in S['K']]); A = np.array([[float(x) for x in r] for r in S['A']]); C = np.array([[float(x) for x in r] for r in S['C']])
L3 = np.array(LAM3, float)
def Evec(ws):
    ws = np.asarray(ws, float); p = np.array([np.sum(ws**m) for m in range(0, 10)])
    E = K + A @ p[1:10]
    pair = np.array([sum(C[i, n - 2] * sum(0.5 * (p[a] * p[n - a] - p[n]) for a in range(1, n)) for n in range(2, 7)) for i in range(8)])
    e3 = (p[1]**3 - 3 * p[1] * p[2] + 2 * p[3]) / 6
    return E + pair + L3 * e3
# check Evec against exact formula on one config
ws = [Fr(1, 3), Fr(2), Fr(5, 7), Fr(1, 9), Fr(3), Fr(7, 2), Fr(1, 2)]
Ex = [S['K'][i] + sum(aval(S, i, x) for x in ws) + sum(bval(S, i, x, y) for x, y in itertools.combinations(ws, 2)) + LAM3[i] * sum(x*y*z for x, y, z in itertools.combinations(ws, 3)) for i in range(8)]
print("Evec check:", np.max(np.abs(Evec([float(x) for x in ws]) - np.array([float(x) for x in Ex])) / (1 + np.abs(np.array([float(x) for x in Ex])))))
rng = np.random.default_rng(7)
cfgs = [np.exp(rng.uniform(-10, 10, 7)) for _ in range(3000)]
cfgs += [np.full(7, np.exp(x)) * np.exp(rng.normal(0, 0.01, 7)) for x in np.linspace(-10, 10, 200)]
for it in range(2):
    R = np.array([Evec(c) for c in cfgs]); N = np.abs(R).sum(axis=1, keepdims=True); keep = (N[:, 0] > 0) & np.all(np.isfinite(R), axis=1); Rn = R[keep] / N[keep]
    # max margin: Rn alpha >= t, |alpha|_1 <= 1
    Aub = np.hstack([-Rn, Rn, np.ones((len(Rn), 1))]); Aub = np.vstack([Aub, np.concatenate([np.ones(16), [0]])])
    bub = np.zeros(len(Aub)); bub[-1] = 1
    r = linprog(np.concatenate([np.zeros(16), [-1]]), A_ub=Aub, b_ub=bub, bounds=[(0, None)] * 16 + [(None, None)], method='highs')
    print('LP status', r.status, r.message, 'rows', len(Rn), 'finite', np.all(np.isfinite(Rn)), flush=True)
    al = r.x[:8] - r.x[8:16]; t = r.x[16]
    # local minimisation of normalised E_alpha from the worst samples
    f = lambda z: (Evec(np.exp(z)) @ al) / np.abs(Evec(np.exp(z))).sum()
    worst = np.argsort(Rn @ al)[:8]; kc = [c for c, kk in zip(cfgs, keep) if kk]; found = []
    for j in worst:
        res = minimize(f, np.log(kc[j]), method='Nelder-Mead', options={'maxiter': 4000, 'xatol': 1e-9, 'fatol': 1e-14})
        found.append((res.fun, np.sort(np.exp(res.x))))
        if np.all(np.abs(res.x) < 30): cfgs.append(np.exp(res.x))
    found.sort(key=lambda x: x[0])
    print(f"iter {it}: margin {t:.3e}  alpha {np.round(al / np.abs(al).max(), 4)}  min found {found[0][0]:.3e} at w = {np.round(found[0][1], 4)}", flush=True)
