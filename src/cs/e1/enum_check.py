"""Independent check of bbc.c (no branch-and-bound, no deletion lemma, no min cut): enumerate EVERY side assignment
sigma in {F, L}^rows, and for each solve the fixed-side closure problem
    max sum x_p - lam sum c_i y_i,  x_p <= y_top(p), y_(q,s,e) <= y_(q,s,e-1), 0 <= x, y <= 1  (sides fixed by sigma)
as an LP with scipy linprog/HiGHS (the closure polytope is integral; integrality of the returned vertex is checked).
M_enum(lam) = max over sigma.  Usage: python3 enum_check.py instfile lam1,lam2,..."""
import sys, itertools
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

def main():
    fn = sys.argv[1]; lams = [float(x) for x in sys.argv[2].split(',')]
    f = open(fn); nP, nI, nR = map(int, f.readline().split())
    items = []
    for _ in range(nI):
        r, s, e, par, c = f.readline().split(); items.append((int(r), int(s), int(e), int(par), float(c)))
    tops = []
    for _ in range(nP):
        a = f.readline().split(); k = int(a[1]); tops.append([int(x) for x in a[2:2 + k]])
    print(f"{fn}: primes {nP}, rows {nR}, items {nI}, side assignments {2 ** nR}", flush=True)
    best = {lam: (-1e300, None) for lam in lams}; nonint = 0
    for sig in itertools.product((0, 1), repeat=nR):
        alive = [i for i in range(nP) if all(items[j][1] == sig[items[j][0]] for j in tops[i])]
        its = [j for j in range(nI) if items[j][1] == sig[items[j][0]]]
        if not alive:
            for lam in lams:
                if 0.0 > best[lam][0]: best[lam] = (0.0, (0, 0.0))
            continue
        ix = {j: k for k, j in enumerate(its)}; nx = len(alive); nv = nx + len(its)
        rr, cc, vv = [], [], []; m = 0
        for k, i in enumerate(alive):
            for j in tops[i]:
                rr += [m, m]; cc += [k, nx + ix[j]]; vv += [1, -1]; m += 1
        for j in its:
            par = items[j][3]
            if par >= 0:
                rr += [m, m]; cc += [nx + ix[j], nx + ix[par]]; vv += [1, -1]; m += 1
        A = coo_matrix((vv, (rr, cc)), shape=(m, nv)).tocsr()
        for lam in lams:
            obj = np.concatenate([-np.ones(nx), np.array([lam * items[j][4] for j in its])])
            res = linprog(obj, A_ub=A, b_ub=np.zeros(m), bounds=(0, 1), method='highs')
            v = -res.fun
            x = res.x
            if np.max(np.minimum(np.abs(x), np.abs(1 - x))) > 1e-7: nonint += 1
            if v > best[lam][0] + 1e-12:
                N = int(round(x[:nx].sum())); b = float(sum(items[j][4] for k, j in enumerate(its) if x[nx + k] > 0.5))
                best[lam] = (v, (N, b))
    for lam in lams:
        v, (N, b) = best[lam]
        print(f"lam {lam:.6f} ENUM M = {v:.6f} at (bits {b:.3f}, N {N})")
    print(f"non-integral LP vertices: {nonint}")

if __name__ == "__main__":
    main()
