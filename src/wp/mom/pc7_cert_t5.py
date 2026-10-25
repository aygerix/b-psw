"""Consistency + diagnosis at shape (1/10, 3/10), D-basis: pair-LP alpha from exact rows; evaluate E_alpha on random
7-member configurations (float, from the symbolic data) and report negatives."""
import sys, itertools, numpy as np
from fractions import Fraction as Fr
from pc7_cert_exact import load as loadx, pair_rows, nrm, lp_feas, grid, aval, bval, LAM3
from pc7_cert_dlib import load as loads, Evec
S = [s for s in loadx('pc7_cert_t3.dat') if s['r'] == (Fr(1, 10), Fr(3, 10))][0]
K, A, C = loads()(0.1, 0.3)
Kx = np.array([float(x) for x in S['K']]); Ax = np.array([[float(x) for x in r] for r in S['A']]); Cx = np.array([[float(x) for x in r] for r in S['C']])
print("symbolic-vs-exact rel diff: K", np.max(np.abs(K - Kx) / (1 + np.abs(Kx))), " A", np.max(np.abs(A - Ax) / (1 + np.abs(Ax))), " C", np.max(np.abs(C - Cx) / (1 + np.abs(Cx))))
G = grid(-14, 14); R = nrm(pair_rows(S, G)); r = lp_feas(R, Kx)
al = r.x; print("pair LP status", r.status, "alpha", np.round(al, 4), "lam3", LAM3 @ al)
rng = np.random.default_rng(0); worst = []
for _ in range(20000):
    ws = np.exp(rng.uniform(-12, 12, 7)); E = Evec(Kx, Ax, Cx, ws); v = E @ al / np.abs(E).sum()
    worst.append((v, ws))
worst.sort(key=lambda x: x[0])
for v, ws in worst[:5]: print(f"  normalised E_alpha {v:+.3e} at w = {np.array2string(np.sort(ws), precision=3)}")
# pair form at the worst config's pairs, and outside the grid range
w0 = np.sort(worst[0][1])
a = lambda x: np.array([np.polyval(np.concatenate([Ax[i][::-1], [0]]), x) for i in range(8)])
b = lambda x, y: np.array([sum(Cx[i, n - 2] * sum(x**q * y**(n - q) for q in range(1, n)) for n in range(2, 7)) for i in range(8)])
pis = [((Kx / 21 + (a(x) + a(y)) / 6 + b(x, y)) @ al, x, y) for x, y in itertools.combinations(w0, 2)]
print("  min pair form over its pairs:", min(pis, key=lambda z: z[0]))
print("  grid range:", float(G[0]), float(G[-1]))
