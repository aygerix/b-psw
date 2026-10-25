#!/usr/bin/env python3
import os
"""src/cs/e3b/advchk/symdp_nodes.py -- node solver on the reflection-symmetric DOUBLE-PENCIL Korselt component (k = 5):
A = x^2 - 1, r_2 = r_3 = r, r_4 = r_5 = -r, l_4 = -l_2, l_5 = -l_3 (f_4(x) = f_2(-x), f_5(x) = f_3(-x)).  This component is
dropped by the E3b nondeg filter (r_i distinct).  Unknowns (c, l_2, l_3, r) + node positions; ONE Korselt equation (K_2;
K_3 = K_2 and K_4, K_5 are mirrors).  Structures list nodes of the independent members 1, 2, 3 (mirrors automatic; member 1's
two roots are +-w and Phi_1 is even, so one node of member 1 counts 2).  Sum m = 6 structures are square; every solution
gets a full Lucas recount (all five members, from the polynomials), so Sum m = 8 points would show up.
Usage: symdp_nodes.py STRUCT nbatch bs seed [--out f]     STRUCT e.g. "1|2|3" (member 1 full, one node each for 2, 3)"""
import sys, argparse, numpy as np
from collections import Counter
sys.path.insert(0, os.path.expanduser('~/b-psw/src/cs/e3b/advchk'))
import sq_enum

class SymDP(sq_enum.Sys):
    def __init__(self, struct, rng=None):
        self.nodes = sq_enum.parse(struct); self.asq = False; self.b = 1.0; self.pencil = None
        self.N = len(self.nodes)
        self.mn = {j: [t for t, S in enumerate(self.nodes) if j in S] for j in range(1, 6)}
        self.nr = 0
        self.n = 4 + self.N
        self.M = self.full(np.full((1, self.n), 0.37 + 0.21j)).shape[1]
        self.R = None
        if self.M > self.n:
            rng = rng or np.random.default_rng(0)
            self.R = rng.normal(size=(self.n, self.M)) + 1j * rng.normal(size=(self.n, self.M))
            self.R[:, :self.n] += 3 * np.eye(self.n)

    def unpack(self, Z):
        c, l2, l3, r = Z[:, 0], Z[:, 1], Z[:, 2], Z[:, 3]
        return c, np.stack([l2, l3, -l2, -l3], axis=1), np.stack([r, r, -r, -r], axis=1), Z[:, 4:]

    def full(self, Z):
        c, l, r, th = self.unpack(Z)
        rj = r[:, 0]; a = rj**2 - 1
        ys = [self.g(c, l, r, i, rj) for i in (1, 3, 4, 5)]
        e = [np.ones_like(c)]
        for y in ys:
            e = [e[0]] + [e[q] + e[q - 1] * y for q in range(1, len(e))] + [e[-1] * y]
        K = np.zeros_like(c)
        for q in range(len(e) - 1, 0, -1):
            K = K * a + e[q]
        res = [K]
        for t, S in enumerate(self.nodes):
            x = th[:, t]
            for j in S:
                pos = self.mn[j].index(t)
                if pos == 0:
                    res.append(self.f(c, l, r, j, x) + 1)
                elif pos == 1:
                    x0 = th[:, self.mn[j][0]]
                    res.append(x0**2 + x0 * x + x**2 - r[:, j - 2] * (x0 + x) - 1)
                else:
                    x0 = th[:, self.mn[j][0]]; x1 = th[:, self.mn[j][1]]
                    res.append(x0 + x1 + x - r[:, j - 2])
            F = np.ones_like(c)
            for i in range(1, 6):
                F = F * self.f(c, l, r, i, x)
            res.append(F + 1)
        return np.stack(res, axis=1)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('struct'); ap.add_argument('nbatch', type=int); ap.add_argument('bs', type=int); ap.add_argument('seed', type=int)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    structs = [l.strip() for l in a.struct.split(',')]
    for st in structs:
        S = SymDP(st, rng=rng)
        fams = {}; conv = 0; hist = []
        for b in range(a.nbatch):
            Z = (rng.normal(size=(a.bs, S.n)) + 1j * rng.normal(size=(a.bs, S.n))) * rng.choice([0.5, 1, 2], size=(a.bs, 1))
            Z, nf = sq_enum.newton(S, Z)
            for s in np.where(nf < 1e-11)[0]:
                z = Z[s]
                if not (np.linalg.norm(S.full(z[None])[0]) < 1e-9):
                    continue
                c, l, r, th = S.unpack(z[None]); c = c[0]; l = l[0]; r = r[0]; th = th[0]
                if abs(c) < 1e-5 or np.min(np.abs(l)) < 1e-5 or abs(l[0] - l[1]) < 1e-5 * max(1, abs(l[0])) or abs(r[0]) < 1e-6 \
                        or np.max(np.abs(z)) > 1e3:
                    continue
                if len(th) > 1 and min(abs(th[u] - th[v]) for u in range(len(th)) for v in range(u + 1, len(th))) < 1e-6:
                    continue
                conv += 1
                kk = sq_enum.famkey(S, z)
                if kk not in fams:
                    tab = sq_enum.lucas(S, z)
                    gd = [sum(1 for (x, rr) in row if rr < 1e-8) for row in tab]
                    fams[kk] = (gd, bool(np.max(np.abs(z[:4].imag)) < 1e-8), z)
            hist.append(len(fams))
        cnt = Counter(tuple(v[0]) for v in fams.values())
        best = max((sum(v[0]) for v in fams.values()), default=None)
        nreal = sum(1 for v in fams.values() if v[1])
        print(f"symdp {st:10s} n={S.n} M={S.M} conv={conv} distinct={len(fams)} real={nreal} max_sum={best} growth={hist[::max(1, len(hist)//6)] + hist[-1:]}  " +
              " ".join(f"{list(k)}x{v}" for k, v in sorted(cnt.items(), key=lambda kv: -sum(kv[0]))[:4]), flush=True)
        if a.out:
            with open(a.out, 'a') as fh:
                for gd, real, z in fams.values():
                    fh.write(repr((st, gd, real, [complex(w) for w in z])) + "\n")

if __name__ == "__main__":
    main()
