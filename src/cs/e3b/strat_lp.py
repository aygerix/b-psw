#!/usr/bin/env python3
"""src/cs/e3b/strat_lp.py -- E3b: Lucas-point parametrisation of the gcd strata (batched LM), any odd k.

A = x^2 - 1, g_1 = c, g_j = l_j (x - r_j).  A member j with m_j >= 1 is parametrised by its first Lucas point th:
  g_j(th) = -2 / A(th)   (i.e. f_j(th) = -1):  c = -2/A(th),  l_j = -2 / (A(th) (th - r_j)).
Members with m_j = 0 carry their own unknown (c or l_j).  Unknowns: r_2..r_k, then per member either its m_j Lucas
points (member 1 with m_1 = 2 uses th and -th, both roots of N_1) or one coefficient.
Equations: Korselt K_j(A(r_j)) = 0 (j >= 2); for every Lucas point th of member j: prod_{i != j} f_i(th) = 1;
for the 2nd/3rd Lucas point of a cubic member: f_j(th) = -1.
Recomputes the true gcd degrees deg gcd(N_j, Phi_j) and the local dimension (rank of the Jacobian) for every
converged non-degenerate solution.
Usage: strat_lp.py k "m1,..,mk" nbatch batchsize seed [--slices d] [--sym] [--out file]
  --sym : restrict to the x -> -x symmetric locus (k = 5: f_3(x) = f_2(-x), f_5(x) = f_4(-x)); m must be symmetric.
"""
import sys, argparse
import numpy as np
from collections import Counter
sys.path.insert(0, __file__.rsplit('/', 1)[0] if '/' in __file__ else '.')
from strat_batch import gcd_degrees, nondeg

class LP:
    def __init__(self, k, m, nslice=0, rng=None, sym=False):
        self.k, self.m, self.sym = k, list(m), sym
        lay = []                       # per member: ('coef', idx) or ('lp', [idx...])
        off = k - 1                    # r's first
        for j in range(k):
            if m[j] == 0:
                lay.append(('coef', off)); off += 1
            elif j == 0:
                lay.append(('lp', [off])); off += 1      # second point is -th
            else:
                lay.append(('lp', list(range(off, off + m[j])))); off += m[j]
        self.lay = lay; self.n = off
        rng = rng or np.random.default_rng(0)
        self.nslice = nslice
        self.S = rng.normal(size=(nslice, self.n)) + 1j * rng.normal(size=(nslice, self.n))
        self.s0 = rng.normal(size=nslice) + 1j * rng.normal(size=nslice)

    def members(self, z):
        k = self.k
        rr = z[:, :k - 1]
        coef = []
        for j in range(k):
            kind, ix = self.lay[j]
            if kind == 'coef':
                coef.append(z[:, ix])
            else:
                th = z[:, ix[0]]; At = th**2 - 1
                coef.append(-2 / At if j == 0 else -2 / (At * (th - rr[:, j - 1])))
        return coef, rr

    def gval(self, coef, rr, j, t):
        return coef[0] if j == 0 else coef[j] * (t - rr[:, j - 1])

    def F(self, z):
        k = self.k
        coef, rr = self.members(z)
        res = []
        for j in range(1, k):
            rj = rr[:, j - 1]; a = rj**2 - 1
            ep = [np.ones_like(rj)]
            for i in range(k):
                if i == j:
                    continue
                y = self.gval(coef, rr, i, rj)
                ep = [ep[0]] + [ep[t] + (ep[t - 1] * y) for t in range(1, len(ep))] + [ep[-1] * y]
            K = np.zeros_like(rj)
            for mm in range(k - 1, 0, -1):
                K = K * a + ep[mm]
            res.append(K)
        for j in range(k):
            kind, ix = self.lay[j]
            if kind != 'lp':
                continue
            pts = [z[:, t] for t in ix]
            if j == 0 and self.m[0] == 2:
                pts = pts + [-pts[0]]
            for s, th in enumerate(pts):
                At = th**2 - 1
                pr = np.ones_like(th)
                for i in range(k):
                    if i != j:
                        pr = pr * (1 + At * self.gval(coef, rr, i, th))
                res.append(pr - 1)
                if s >= 1 and j >= 1:
                    res.append(2 + At * self.gval(coef, rr, j, th))
        if self.sym:          # k = 5: member 3 = mirror of 2, member 5 = mirror of 4: l3 = -l2, r3 = -r2, ...
            res.append(rr[:, 1] + rr[:, 0]); res.append(coef[2] + coef[1])
            res.append(rr[:, 3] + rr[:, 2]); res.append(coef[4] + coef[3])
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

def lm_batch(model, z, iters=200, tol=1e-13, prune=False):
    B, n = z.shape
    with np.errstate(all='ignore'):
        f = model.F(z)
    nf = np.linalg.norm(f, axis=1)
    active = np.isfinite(nf)
    z[~active] = 0.3; f[~active] = 0; nf[~active] = np.inf
    mu = np.full(B, 1e-3)
    def hess(idx):
        with np.errstate(all='ignore'):
            Jz = model.J(z[idx])
        Jz[~np.isfinite(Jz)] = 0
        JH = np.conj(np.transpose(Jz, (0, 2, 1)))
        return JH @ Jz, (JH @ f[idx][:, :, None])[:, :, 0]
    H = np.zeros((B, n, n), dtype=complex); g = np.zeros((B, n), dtype=complex)
    i0 = np.where(active)[0]
    if len(i0):
        H[i0], g[i0] = hess(i0)
    eye = np.eye(n)[None]
    for it in range(iters):
        idx = np.where(active & (nf > tol))[0]
        if len(idx) == 0:
            break
        try:
            d = np.linalg.solve(H[idx] + mu[idx][:, None, None] * eye, -g[idx][:, :, None])[:, :, 0]
        except np.linalg.LinAlgError:
            d = np.zeros((len(idx), n), dtype=complex)
        zn = z[idx] + d
        with np.errstate(all='ignore'):
            fn = model.F(zn); nn = np.linalg.norm(fn, axis=1)
        ok = np.isfinite(nn) & (nn < nf[idx])
        acc = idx[ok]; rej = idx[~ok]
        z[acc] = zn[ok]; f[acc] = fn[ok]; nf[acc] = nn[ok]
        mu[acc] = np.maximum(mu[acc] / 3, 1e-12); mu[rej] *= 4
        if len(acc):
            H[acc], g[acc] = hess(acc)
        active &= (mu < 1e10) & (np.max(np.abs(z), axis=1) < 1e6)
        if prune and it == 60:
            active &= nf < 1e-3
    return z, nf

def to_clr(model, zs):
    coef, rr = model.members(zs[None])
    return np.array([c[0] for c in coef] + list(rr[0]))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('k', type=int); ap.add_argument('m'); ap.add_argument('nbatch', type=int)
    ap.add_argument('bs', type=int); ap.add_argument('seed', type=int)
    ap.add_argument('--slices', type=int, default=0); ap.add_argument('--out', default=None)
    ap.add_argument('--sym', action='store_true'); ap.add_argument('--iters', type=int, default=200)
    a = ap.parse_args()
    k = a.k; m = [int(t) for t in a.m.split(',')]
    rng = np.random.default_rng(a.seed)
    model = LP(k, m, nslice=a.slices, rng=rng, sym=a.sym)
    M = model.F(np.full((1, model.n), 0.3 + 0.1j)).shape[1]
    print(f"k={k} m={m} sym={a.sym} slices={a.slices} unknowns={model.n} equations={M} expected dim={model.n - M}", flush=True)
    found = []
    for b in range(a.nbatch):
        scale = rng.choice([0.5, 1.0, 2.0, 4.0], size=(a.bs, 1))
        z = (rng.normal(size=(a.bs, model.n)) + 1j * rng.normal(size=(a.bs, model.n))) * scale
        z, nf = lm_batch(model, z, iters=a.iters)
        for s in np.where(nf < 1e-11)[0]:
            zs = z[s]
            clr = to_clr(model, zs)
            if not np.all(np.isfinite(clr)) or not nondeg(clr, k):
                continue
            # distinct Lucas points within a member
            okd = True
            for j in range(k):
                kind, ix = model.lay[j]
                if kind == 'lp' and len(ix) > 1:
                    pts = zs[ix]
                    for u in range(len(pts)):
                        for v in range(u + 1, len(pts)):
                            if abs(pts[u] - pts[v]) < 1e-6:
                                okd = False
            if not okd:
                continue
            Jz = model.J(zs[None])[0]
            sv = np.linalg.svd(Jz, compute_uv=False)
            rank = int(np.sum(sv > 1e-7 * sv[0]))
            gd, pts = gcd_degrees(clr, k)
            found.append((model.n - rank, gd, clr, pts))
        print(f"  batch {b}: kept {len(found)}", flush=True)
    cnt = Counter((d, tuple(gd)) for d, gd, clr, pts in found)
    print(f"converged non-degenerate: {len(found)} / {a.nbatch * a.bs}")
    for key, v in sorted(cnt.items()):
        print(f"  local dim {key[0]}  actual gcd degrees {list(key[1])} (sum {sum(key[1])}): {v}")
    if a.out:
        with open(a.out, 'w') as fh:
            for d, gd, clr, pts in found:
                fh.write(repr((d, gd, [complex(t) for t in clr], [[complex(t) for t in p] for p in pts])) + "\n")

if __name__ == "__main__":
    main()
