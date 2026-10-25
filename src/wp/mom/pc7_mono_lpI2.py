"""Termwise LP on the reality interval, row-normalised (scale-free margin); reports binding rows."""
import numpy as np, sys
from scipy.optimize import linprog
from pc7_mono_lib import parts
from pc7_mono_lpI import interval
def build(r1, r2, ns=400):
    I = interval(r1, r2)
    if I is None: return None
    lo, hi = I; hi_s = min(hi, lo * 1e4)
    X = np.exp(np.linspace(np.log(lo), np.log(hi_s), ns))
    C0, P1, A, kap = map(np.array, parts(r1, r2))
    R1 = np.array([[sum(P1[i][a] * x**(a - 1) for a in range(1, 10)) for i in range(8)] for x in X])
    R2 = np.array([[sum((m + 1) * A[i][m] * x**m for m in range(5)) for i in range(8)] for x in X])
    return I, X, C0, R1, R2, kap
def lp(r1, r2, supp=range(8), drop=()):
    b = build(r1, r2)
    if b is None: return 'EMPTY', None, None, None
    I, X, C0, R1, R2, kap = b; supp = list(supp); n = len(supp)
    blocks = {'phi1': R1, 'gam': R2}
    rows, tags = [], []
    for nm, R in blocks.items():
        if nm in drop: continue
        for s, x in enumerate(X):
            v = R[s, supp]; nv = np.linalg.norm(v)
            if nv > 0: rows.append(v / nv); tags.append((nm, x))
    M = np.array(rows)
    Aub = np.hstack([-M, np.ones((len(M), 1))]); bub = np.zeros(len(M))
    kv = kap[supp]
    if np.any(kv != 0): Aub = np.vstack([Aub, np.concatenate([-kv, [0]])]); bub = np.append(bub, 0)
    c0 = C0[supp]; Aeq = np.concatenate([c0 / np.linalg.norm(c0), [0]])[None, :]
    r = linprog(np.concatenate([np.zeros(n), [-1]]), A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=[1],
                bounds=[(-1e3, 1e3)] * n + [(None, 1)], method='highs')
    if r.status != 0: return 'LPFAIL', I, None, None
    al = r.x[:n]; sl = M @ al
    bind = [tags[i] for i in np.argsort(sl)[:4]]
    return r.x[-1], I, al, bind
if __name__ == '__main__':
    np.set_printoptions(linewidth=200, precision=4)
    for r1, r2 in [(-0.3, 0.6), (0.1, 0.5), (0.0, 0.2), (-0.5, 0.5), (-0.1, 0.1), (-0.9, 0.2), (0.3, 0.9), (-0.95, 0.9), (0.0, 0.05)]:
        d, I, a, bind = lp(r1, r2)
        print(f"r = ({r1:+.2f},{r2:+.2f}) I = {None if I is None else np.round(I, 5)} margin {d}\n   alpha {a}\n   binding {bind}", flush=True)
        for drop in [('gam',), ('phi1',)]:
            d2, _, a2, b2 = lp(r1, r2, drop=drop)
            print(f"   without {drop}: margin {d2}")
