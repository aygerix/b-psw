"""Full 7-member LP in the D-basis with robust sampling, at given shapes, for several alpha-subspaces:
  all : alpha in R^8;   S : alpha_i = alpha_(i+4) (sums D_i + D_(i+4): top power sums and e3 cancel);  V, W : one side only.
Rows: random log-uniform configurations w in [1e-6, 1e6]^7, clustered configurations (1..3 distinct values), plus the
w_j -> infinity asymptotic rows (coefficients of w_j^9, w_j^8, w_j^7, w_j^6 are shape constants, no dependence on the others).
Prints max margin t (rows normalised, |alpha|_1 <= 1) and the optimal alpha.  Usage: python3 pc7_cert_dlp.py nconf r1,r2 [r1,r2 ...]"""
import sys, numpy as np
from scipy.optimize import linprog
from pc7_cert_dlib import load, Evec, LAM3
at = load(); rng = np.random.default_rng(5); nconf = int(sys.argv[1])
def configs():
    cf = [np.exp(rng.uniform(-14, 14, 7)) for _ in range(nconf)]
    for _ in range(nconf // 2):
        m = rng.integers(1, 4); vals = np.exp(rng.uniform(-14, 14, m)); idx = rng.integers(0, m, 7)
        cf.append(vals[idx] * np.exp(rng.normal(0, 1e-3, 7)))
    return cf
CF = configs()
SUB = {'all': np.eye(8), 'D1': np.eye(8)[:, :1], 'D15': np.eye(8)[:, [0, 4]], 'S': np.vstack([np.eye(4), np.eye(4)]), 'V': np.vstack([np.eye(4), np.zeros((4, 4))]), 'W': np.vstack([np.zeros((4, 4)), np.eye(4)])}
for arg in sys.argv[2:]:
    r1, r2 = map(float, arg.split(','))
    K, A, C = at(r1, r2); E1 = r1 + r2; lam = 2 * E1 / ((1 - r1**2) * (1 - r2**2))
    R = np.array([Evec(K, A, C, c) for c in CF]); R = R / np.abs(R).sum(axis=1, keepdims=True)
    for name, M in SUB.items():
        Rm = R @ M; d = M.shape[1]
        # asymptotic: the w^9 row (lam^8 * LAM3) must be >= 0; if it vanishes for alpha, w^8 must be >= 0 etc. (lexicographic):
        # here impose lam3 >= 0 and A_alpha,8 >= 0 when lam3 = 0 is forced (S); a sufficient linear version: require each >= 0
        asym = [LAM3 @ M] + ([A[:, 7] @ M] if name == 'S' else [])
        Aub = np.vstack([np.hstack([-Rm, Rm, np.ones((len(Rm), 1))])] + [np.concatenate([-a, a, [0]])[None, :] for a in asym] + [np.concatenate([np.ones(2 * d), [0]])[None, :]])
        bub = np.zeros(len(Aub)); bub[-1] = 1
        r = linprog(np.concatenate([np.zeros(2 * d), [-1]]), A_ub=Aub, b_ub=bub, bounds=[(0, None)] * (2 * d) + [(None, None)], method='highs')
        if r.status != 0: print(f"r = ({r1:+.3f},{r2:+.3f}) {name:3s}: LP status {r.status}"); continue
        al = M @ (r.x[:d] - r.x[d:2 * d]); t = r.x[-1]
        print(f"r = ({r1:+.3f},{r2:+.3f}) {name:3s}: margin {t:+.3e}  alpha {np.array2string(al / max(np.abs(al).max(), 1e-300), precision=4, suppress_small=True)}", flush=True)
