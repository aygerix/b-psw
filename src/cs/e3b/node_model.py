#!/usr/bin/env python3
"""src/cs/e3b/node_model.py -- E3b: Lucas NODE model with explicit sharing structure (k = 5), batched LM.

A = x^2 - 1, g_1 = c, g_j = l_j (x - r_j).  An incidence structure is a list of nodes; node t has a member set S_t
(f_j(th_t) = -1 for j in S_t).  sum |S_t| = Sum m.  Member j with >= 1 node is parametrised by its FIRST node th:
  member 1: c = -2/A(th);  cubic j: f_j(x) = 1 - 2 E_j(x)/E_j(th), E_j(x) = (x^2 - 1)(x - r_j).
Equations: Korselt_j (j = 2..5);  further nodes of member j: member 1: th' = -th;  cubic: divided differences
[E_j](th1, th2) = 0, [E_j](th1, th2, th3) = th1 + th2 + th3 - r_j = 0  (these exclude coincidences);
per node with |S| < 5: prod_{i not in S} f_i(th) = (-1)^(|S| - 1).
Usage: node_model.py STRUCT nbatch bs seed [--out f]   STRUCT like "1,1|2|3|45" ('|' separates nodes; digits = members)
       node_model.py --enum SUM   (print canonical structures with Sum m = SUM up to permutations of members 2..5)
"""
import sys, argparse, itertools
import numpy as np
from collections import Counter
sys.path.insert(0, __file__.rsplit('/', 1)[0] if '/' in __file__ else '.')
from strat_lp import lm_batch
from refine import lucas_table
from strat_batch import nondeg

CAP = {1: 2, 2: 3, 3: 3, 4: 3, 5: 3}

def parse(s):
    return [tuple(sorted(int(ch) for ch in node if ch.isdigit())) for node in s.split('|')]

def canon(nodes):
    best = None
    for perm in itertools.permutations([2, 3, 4, 5]):
        mp = {1: 1, 2: perm[0], 3: perm[1], 4: perm[2], 5: perm[3]}
        nn = sorted(tuple(sorted(mp[j] for j in node)) for node in nodes)
        key = tuple(nn)
        if best is None or key < best:
            best = key
    return best

def enum(total):
    subsets = [tuple(c) for r in range(1, 6) for c in itertools.combinations(range(1, 6), r)]
    out = set()
    def rec(start, nodes, cnt, s):
        if s == total:
            out.add(canon(nodes)); return
        for i in range(start, len(subsets)):
            S = subsets[i]
            if s + len(S) > total:
                continue
            if any(cnt[j] + 1 > CAP[j] for j in S):
                continue
            for j in S: cnt[j] += 1
            rec(i, nodes + [S], cnt, s + len(S))
            for j in S: cnt[j] -= 1
    rec(0, [], {j: 0 for j in range(1, 6)}, 0)
    return sorted(out)

class Node:
    def __init__(self, nodes):
        self.nodes = nodes
        self.first = {}; self.mnodes = {j: [] for j in range(1, 6)}
        for t, S in enumerate(nodes):
            for j in S:
                self.mnodes[j].append(t)
        # unknown layout: r2..r5 (0..3), node positions (4..4+N-1), coefficients for node-less members
        self.N = len(nodes)
        off = 4 + self.N
        self.coef = {}
        for j in range(1, 6):
            if not self.mnodes[j]:
                self.coef[j] = off; off += 1
        self.n = off

    def fval(self, z, j, x):
        """f_j(x) for batch z."""
        if self.mnodes[j]:
            th = z[:, 4 + self.mnodes[j][0]]
            if j == 1:
                return 1 - 2 * (x**2 - 1) / (th**2 - 1)
            r = z[:, j - 2]
            return 1 - 2 * (x**2 - 1) * (x - r) / ((th**2 - 1) * (th - r))
        cc = z[:, self.coef[j]]
        if j == 1:
            return 1 + (x**2 - 1) * cc
        return 1 + (x**2 - 1) * cc * (x - z[:, j - 2])

    def gval(self, z, j, x):
        if self.mnodes[j]:
            th = z[:, 4 + self.mnodes[j][0]]
            if j == 1:
                return -2 / (th**2 - 1) + 0 * x
            r = z[:, j - 2]
            return -2 * (x - r) / ((th**2 - 1) * (th - r))
        cc = z[:, self.coef[j]]
        return cc + 0 * x if j == 1 else cc * (x - z[:, j - 2])

    def F(self, z):
        res = []
        for j in range(2, 6):
            rj = z[:, j - 2]; a = rj**2 - 1
            ep = [np.ones_like(rj)]
            for i in range(1, 6):
                if i == j:
                    continue
                y = self.gval(z, i, rj)
                ep = [ep[0]] + [ep[q] + ep[q - 1] * y for q in range(1, len(ep))] + [ep[-1] * y]
            K = np.zeros_like(rj)
            for mm in range(4, 0, -1):
                K = K * a + ep[mm]
            res.append(K)
        for j in range(1, 6):
            ts = self.mnodes[j]
            if len(ts) <= 1:
                continue
            th1 = z[:, 4 + ts[0]]
            if j == 1:
                res.append(z[:, 4 + ts[1]] + th1)
                continue
            r = z[:, j - 2]
            for u in ts[1:]:
                th = z[:, 4 + u]
                res.append(th**2 + th * th1 + th1**2 - 1 - r * (th + th1))
            if len(ts) == 3:
                res.append(th1 + z[:, 4 + ts[1]] + z[:, 4 + ts[2]] - r)
        for t, S in enumerate(self.nodes):
            if len(S) == 5:
                continue
            th = z[:, 4 + t]
            pr = np.ones_like(th)
            for i in range(1, 6):
                if i not in S:
                    pr = pr * self.fval(z, i, th)
            res.append(pr - (-1)**(len(S) - 1))
        return np.stack(res, axis=1)

    def J(self, z, h=1e-6):
        B, n = z.shape
        cols = []
        for i in range(n):
            e = np.zeros(n); e[i] = 1
            hh = h * np.maximum(1, np.abs(z[:, i]))[:, None]
            cols.append((self.F(z + hh * e) - self.F(z - hh * e)) / (2 * hh))
        return np.stack(cols, axis=2)

    def clr(self, zs):
        z = zs[None]
        c = self.gval(z, 1, 0.0)[0]
        lam = []
        for j in range(2, 6):
            g0 = self.gval(z, j, 0.0)[0]; g1 = self.gval(z, j, 1.0)[0]
            lam.append(g1 - g0)
        return np.array([c] + lam + list(zs[:4]))

def run(struct, nbatch, bs, seed, out=None, quiet=False, minsum=7):
    nodes = parse(struct)
    model = Node(nodes)
    rng = np.random.default_rng(seed)
    M = model.F(np.full((1, model.n), 0.3 + 0.2j)).shape[1]
    found = []
    for b in range(nbatch):
        z = (rng.normal(size=(bs, model.n)) + 1j * rng.normal(size=(bs, model.n))) * rng.choice([0.5, 1, 2], size=(bs, 1))
        with np.errstate(all='ignore'):
            z, nf = lm_batch(model, z, iters=200, prune=True)
        for s in np.where(nf < 1e-11)[0]:
            zs = z[s]
            th = zs[4:4 + model.N]
            if model.N > 1 and min(abs(th[u] - th[v]) for u in range(model.N) for v in range(u + 1, model.N)) < 1e-6:
                continue                        # nodes must be distinct points
            clr = model.clr(zs)
            if not np.all(np.isfinite(clr)) or not nondeg(clr, 5):
                continue
            tab = lucas_table(clr, 5)
            gd = [sum(1 for (t, r) in row if r < 1e-7) for row in tab]
            found.append((gd, clr, zs))
    cnt = Counter(tuple(gd) for gd, clr, zs in found)
    best = max([sum(g) for g in cnt], default=None)
    line = f"{struct:24s} n={model.n} M={M} exp={model.n - M:3d} conv={len(found):5d}/{nbatch * bs}  max sum m = {best}  " + \
           " ".join(f"{list(k)}x{v}" for k, v in sorted(cnt.items(), key=lambda kv: -sum(kv[0]))[:3])
    print(line, flush=True)
    if out:
        with open(out, 'a') as fh:
            for gd, clr, zs in found:
                if sum(gd) >= minsum:
                    fh.write(repr((struct, gd, [complex(u) for u in clr], [complex(u) for u in zs])) + "\n")
    return best

def main():
    if sys.argv[1] == '--enum':
        for s in enum(int(sys.argv[2])):
            print('|'.join(''.join(str(j) for j in S) for S in s))
        return
    ap = argparse.ArgumentParser()
    ap.add_argument('struct'); ap.add_argument('nbatch', type=int); ap.add_argument('bs', type=int); ap.add_argument('seed', type=int)
    ap.add_argument('--out', default=None); ap.add_argument('--minsum', type=int, default=7)
    a = ap.parse_args()
    if a.struct.startswith('@'):
        for line in open(a.struct[1:]):
            if line.strip():
                run(line.strip(), a.nbatch, a.bs, a.seed, a.out, minsum=a.minsum)
    else:
        run(a.struct, a.nbatch, a.bs, a.seed, a.out, minsum=a.minsum)

if __name__ == "__main__":
    main()
