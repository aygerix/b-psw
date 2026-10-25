"""Structure: b_i(w, v) = w v sum_{N=0}^{4} c_{i,N} h_N(w, v) (complete homogeneous), so the pair part is a Hankel form.
Tests at grid shapes (all sufficient for sum alpha_i E_i > 0 on (0, inf)^7):
  dp7 : kap = 1, c_alpha,N >= 0 (N = 0..4) [so b_alpha >= 0], a_alpha(w) >= -1/7 (w > 0), lam3 >= 0
  b0  : kap = 1, c_alpha,N = 0 [b_alpha == 0], a_alpha(w) >= -1/7, lam3 >= 0
Usage: python3 pc7_cert_dp7.py ngrid"""
import sys, numpy as np
from scipy.optimize import linprog
from pc7_cert_lib import load, at_shape, LAM3, NAMES
D = load()
WG = np.exp(np.linspace(np.log(1e-6), np.log(1e6), 301))
def hstruct(Bc):
    """check b_i[m][n] depends only on m + n; return c[i, N] = coeff of w v^(N+1)"""
    err = 0.0; C = np.zeros((8, 5))
    for i in range(8):
        for N in range(5): C[i, N] = Bc[i, 1, N + 1]
        for m in range(10):
            for n in range(10):
                tgt = C[i, m + n - 2] if (m >= 1 and n >= 1 and m + n - 2 <= 4) else 0.0
                err = max(err, abs(Bc[i, m, n] - tgt) / (1 + abs(tgt)))
    return C, err
def test(R1, R2, mode):
    kap, Ac, Bc = at_shape(D, R1, R2); C, err = hstruct(Bc)
    Ra = np.array([[np.polyval(Ac[i][::-1], x) for i in range(8)] for x in WG]); sc = np.abs(Ra).max(axis=1)
    # a_alpha(w) + kap_alpha/7 >= 0, scaled
    Aub = [-(Ra + kap[None, :] / 7) / np.maximum(sc, 1)[:, None], -LAM3[None, :]]
    Aeq = [kap[None, :]]; beq = [1.0]
    if mode == 'dp7': Aub.append(-C.T / np.maximum(np.abs(C.T).max(axis=1, keepdims=True), 1e-300))
    else: Aeq.append(C.T); beq += [0.0] * 5
    r = linprog(np.zeros(8), A_ub=np.vstack(Aub), b_ub=np.zeros(sum(len(x) for x in Aub)), A_eq=np.vstack(Aeq), b_eq=beq,
                bounds=[(-1e5, 1e5)] * 8, method='highs')
    return r.status == 0, err
ng = int(sys.argv[1]); grid = np.linspace(-0.9, 0.9, ng)
shapes = [(a, b) for a in grid for b in grid if b >= a - 1e-12]
for mode in ('dp7', 'b0'):
    res = [(s, *test(*s, mode)) for s in shapes]
    bad = [(round(s[0], 2), round(s[1], 2)) for s, ok, e in res if not ok]
    print(f"mode {mode}: {len(res) - len(bad)}/{len(res)} shapes certified; max Hankel-structure error {max(e for _, _, e in res):.1e}; failing {bad}", flush=True)
