#!/usr/bin/env python3
"""src/cs/e3b/fsp_num.py -- E3b: numerical search INSIDE the FSP locus (k = 5), clean coordinates.
A = x^2 - 1, E_i(x) = (x^2 - 1)(x - r_i); f_1 = 1 - 2 (x^2 - 1)/(t^2 - 1), f_i = 1 - 2 E_i(x)/E_i(t)  (i = 2..5).
Unknowns t, r_2..r_5 (+ extra Lucas points).  Korselt_j: prod_{i != j} f_i(r_j) = 1 (K-form, divided by r_j^2 - 1).
Extra incidences (list): '1' = f_1 full (Lucas at -t); 'j' (j = 2..5) = an extra Lucas point s of member j:
   Q_j(s) = s^2 + (t - r_j) s + t^2 - r_j t - 1 = 0 and prod_{i != j} f_i(s) = 1;  'jj' = both extra points of j.
Usage: fsp_num.py "1,2" nbatch bs seed [--slices d]
Reports converged non-degenerate solutions with their recomputed gcd degrees (refine.lucas_table) and local dimension.
"""
import sys, argparse
import numpy as np
from collections import Counter
sys.path.insert(0, __file__.rsplit('/', 1)[0] if '/' in __file__ else '.')
from strat_lp import lm_batch
from refine import lucas_table

class FSP:
    def __init__(self, extra, nslice=0, rng=None):
        self.extra = extra
        self.n = 5 + sum(1 if e != '1' and len(e) == 1 else (2 if len(e) == 2 else 0) for e in extra)
        rng = rng or np.random.default_rng(0)
        self.nslice = nslice
        self.S = rng.normal(size=(nslice, self.n)) + 1j * rng.normal(size=(nslice, self.n))
        self.s0 = rng.normal(size=nslice) + 1j * rng.normal(size=nslice)

    def f(self, z, i, x):
        t = z[:, 0]
        if i == 1:
            return 1 - 2 * (x**2 - 1) / (t**2 - 1)
        r = z[:, i - 1]
        return 1 - 2 * (x**2 - 1) * (x - r) / ((t**2 - 1) * (t - r))

    def F(self, z):
        t = z[:, 0]
        res = []
        for j in range(2, 6):
            rj = z[:, j - 1]; a = rj**2 - 1
            ys = []
            for i in range(1, 6):
                if i == j:
                    continue
                if i == 1:
                    ys.append(-2 / (t**2 - 1))
                else:
                    ri = z[:, i - 1]
                    ys.append(-2 * (rj - ri) / ((t**2 - 1) * (t - ri)))
            ep = [np.ones_like(t)]
            for y in ys:
                ep = [ep[0]] + [ep[q] + ep[q - 1] * y for q in range(1, len(ep))] + [ep[-1] * y]
            K = np.zeros_like(t)
            for mm in range(4, 0, -1):
                K = K * a + ep[mm]
            res.append(K)
        off = 5
        for e in self.extra:
            if e == '1':
                pr = np.ones_like(t)
                for i in range(2, 6):
                    pr = pr * self.f(z, i, -t)
                res.append(pr - 1)
            else:
                j = int(e[0])
                rj = z[:, j - 1]
                for s_ in range(len(e)):
                    s = z[:, off]; off += 1
                    res.append(s**2 + (t - rj) * s + t**2 - rj * t - 1)
                    pr = np.ones_like(t)
                    for i in range(1, 6):
                        if i != j:
                            pr = pr * self.f(z, i, s)
                    res.append(pr - 1)
                if len(e) == 2:   # the two extra points must be the two distinct roots: sum = r_j - t
                    res.append(z[:, off - 2] + z[:, off - 1] - (rj - t))
        out = np.stack(res, axis=1)
        if self.nslice:
            out = np.concatenate([out, z @ self.S.T - self.s0], axis=1)
        return out

    def J(self, z, h=1e-6):
        B, n = z.shape
        cols = []
        for i in range(n):
            e = np.zeros(n); e[i] = 1
            hh = h * np.maximum(1, np.abs(z[:, i]))[:, None]
            cols.append((self.F(z + hh * e) - self.F(z - hh * e)) / (2 * hh))
        return np.stack(cols, axis=2)

def to_clr(z):
    t = z[0]; rr = z[1:5]
    c = -2 / (t**2 - 1)
    lam = [c / (t - r) for r in rr]
    return np.array([c] + lam + list(rr))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('extra'); ap.add_argument('nbatch', type=int); ap.add_argument('bs', type=int)
    ap.add_argument('seed', type=int); ap.add_argument('--slices', type=int, default=0); ap.add_argument('--out', default=None)
    a = ap.parse_args()
    extra = [e for e in a.extra.split(',') if e]
    rng = np.random.default_rng(a.seed)
    model = FSP(extra, nslice=a.slices, rng=rng)
    M = model.F(np.full((1, model.n), 0.3 + 0.2j)).shape[1]
    print(f"FSP extra={extra} slices={a.slices} unknowns={model.n} equations={M} expected dim={model.n - M}", flush=True)
    found = []
    for b in range(a.nbatch):
        z = (rng.normal(size=(a.bs, model.n)) + 1j * rng.normal(size=(a.bs, model.n))) * rng.choice([0.5, 1, 2], size=(a.bs, 1))
        z, nf = lm_batch(model, z, iters=200)
        for s in np.where(nf < 1e-11)[0]:
            zs = z[s]; t = zs[0]; rr = zs[1:5]
            if abs(t**2 - 1) < 1e-5 or min(abs(t - r) for r in rr) < 1e-5 or np.max(np.abs(zs)) > 300:
                continue
            if min(abs(rr[i] - rr[j]) for i in range(4) for j in range(i + 1, 4)) < 1e-5:
                continue
            clr = to_clr(zs)
            tab = lucas_table(clr, 5)
            gd = [sum(1 for (tt, r) in row if r < 1e-7) for row in tab]
            Jz = model.J(zs[None])[0]
            sv = np.linalg.svd(Jz, compute_uv=False)
            rank = int(np.sum(sv > 1e-7 * sv[0]))
            found.append((model.n - rank, gd, clr, zs))
    cnt = Counter((d, tuple(gd)) for d, gd, clr, zs in found)
    print(f"converged non-degenerate: {len(found)} / {a.nbatch * a.bs}")
    for key, v in sorted(cnt.items()):
        print(f"  local dim {key[0]}  gcd degrees {list(key[1])} (sum {sum(key[1])}): {v}")
    distinct = {}
    for d, gd, clr, zs in found:
        key = tuple(np.round(sorted([abs(zs[0])] + list(np.abs(zs[1:5]))), 5))
        distinct.setdefault(key, (d, gd, zs))
    print("distinct points (up to x -> -x):", len(distinct))
    if a.out:
        with open(a.out, 'w') as fh:
            for key, (d, gd, zs) in distinct.items():
                fh.write(repr((d, gd, [complex(u) for u in zs])) + "\n")

if __name__ == "__main__":
    main()
