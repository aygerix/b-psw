"""E1 direct ILP (cs_rank.md formulation) with scipy.optimize.milp (HiGHS), per (D, class, b):
  max sum_p x_p   s.t.  x_p <= y_{q,s,a} for every requirement (q, s, a) of p;  y_{q,s,e} <= y_{q,s,e-1};
                        y_{q,F,1} + y_{q,L,1} <= 1;  sum_{q,s,e} c(q,e) y_{q,s,e} <= b;  x, y binary.
Reports the LP-relaxation bound, the MILP dual bound at the time limit (both valid upper bounds), and the best
integer solution.  Usage: python3 milp_e1.py D j1 j2 b1,b2,... TLIM recfile1 [recfile2 ...]
"""
import sys, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix
sys.path.insert(0, '.')
from lagr_bb import load, Inst, cost

def build(inst, b):
    nP, nI = inst.n, len(inst.items)
    nv = nP + nI
    r, c, v = [], [], []; lo, hi = [], []; m = 0
    for i, req in enumerate(inst.req):
        for q, s, a in req:
            r += [m, m]; c += [i, nP + inst.items[(q, s, a)]]; v += [1, -1]; lo.append(-np.inf); hi.append(0); m += 1
    for (q, s, e), j in inst.items.items():
        if e > 1:
            r += [m, m]; c += [nP + j, nP + inst.items[(q, s, e - 1)]]; v += [1, -1]; lo.append(-np.inf); hi.append(0); m += 1
    for q in inst.rows:
        ks = [inst.items[(q, s, 1)] for s in (0, 1) if (q, s, 1) in inst.items]
        if len(ks) == 2:
            r += [m, m]; c += [nP + ks[0], nP + ks[1]]; v += [1, 1]; lo.append(-np.inf); hi.append(1); m += 1
    for (q, s, e), j in inst.items.items():
        r.append(m); c.append(nP + j); v.append(cost(q, e))
    lo.append(-np.inf); hi.append(b); m += 1
    A = coo_matrix((v, (r, c)), shape=(m, nv)).tocsr()
    obj = np.concatenate([-np.ones(nP), np.zeros(nI)])
    return obj, LinearConstraint(A, lo, hi), Bounds(0, 1)

def main():
    D, j1, j2 = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    bl = [float(x) for x in sys.argv[4].split(',')]; tl = float(sys.argv[5])
    pool, bad = load(D, j1, j2, sys.argv[6:])
    inst = Inst(pool)
    print(f"D = {D}, class ({j1},{j2}): pool {inst.n}, items {len(inst.items)}, rows {len(inst.rows)}", flush=True)
    for b in bl:
        obj, con, bnd = build(inst, b)
        t0 = time.time()
        lp = milp(obj, constraints=con, bounds=bnd, integrality=np.zeros(len(obj)))
        t1 = time.time()
        ip = milp(obj, constraints=con, bounds=bnd, integrality=np.ones(len(obj)),
                  options=dict(time_limit=tl, disp=False, mip_rel_gap=0))
        t2 = time.time()
        best = -ip.fun if ip.x is not None else None
        dual = -ip.mip_dual_bound if getattr(ip, 'mip_dual_bound', None) is not None else None
        print(f"  b = {b:6.1f}: LP bound {-lp.fun:9.2f} ({t1 - t0:.1f} s) | MILP status {ip.status} "
              f"best integer {best} dual bound {dual if dual is None else round(dual, 2)} "
              f"nodes {getattr(ip, 'mip_node_count', None)} ({t2 - t1:.0f} s) | (b+1)/3 = {(b + 1) / 3:.1f}", flush=True)

if __name__ == "__main__":
    main()
