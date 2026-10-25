#!/usr/bin/env python3
"""src/cs/e3b/strat_batch.py -- E3b: batched numerical survey of Lucas gcd-degree strata (any odd k).

Same model as strat_num.py (A = x^2 - 1, g_1 = c, g_j = l_j (x - r_j)), evaluated with vectorised numpy over a batch
of random starts (complex Levenberg-Marquardt, central-difference Jacobian).  Polynomials are coefficient arrays,
lowest degree first, shape (B, deg + 1).

Pattern m = (m_1..m_k): deg gcd(N_j, Phi_j) >= m_j, via monic Q_j (deg m_j) with N_j == Phi_j == 0 mod Q_j,
or Phi_j == 0 mod N_j when m_j = deg N_j.
Optional extra constraints (flags):
  --fsp     : a fully shared Lucas point t (unknown): f_j(t) = -1 for all j (equivalently g_j(t) = c for all j and
              c A(t) = -2); this is only used to steer the search, the reported gcd degrees are recomputed from scratch.
  --slices d: d random affine slices (to sample a d-dimensional component).
Usage: strat_batch.py k "m1,..,mk" nbatch batchsize seed [--fsp] [--slices d] [--out file]
"""
import sys, argparse
import numpy as np
from collections import Counter

def pmul(a, b):
    B = a.shape[0]
    out = np.zeros((B, a.shape[1] + b.shape[1] - 1), dtype=complex)
    for i in range(a.shape[1]):
        out[:, i:i + b.shape[1]] += a[:, i:i + 1] * b
    return out

def padd(a, b):
    n = max(a.shape[1], b.shape[1])
    out = np.zeros((a.shape[0], n), dtype=complex)
    out[:, :a.shape[1]] += a
    out[:, :b.shape[1]] += b
    return out

def prem_monic(P, Q):
    """remainder of P by monic Q (Q given as (B, m+1) with Q[:, m] == 1)."""
    m = Q.shape[1] - 1
    R = P.copy()
    for t in range(R.shape[1] - 1, m - 1, -1):
        co = R[:, t:t + 1]
        R[:, t - m:t + 1] -= co * Q
    return R[:, :m]

class Model:
    def __init__(self, k, m, fsp=False, nslice=0, rng=None):
        self.k, self.m, self.fsp, self.nslice = k, list(m), fsp, nslice
        self.nbase = 2 * k - 1
        self.qoff = []
        off = self.nbase
        for j in range(k):
            dN = 2 if j == 0 else 3
            if 0 < m[j] < dN:
                self.qoff.append((j, off, m[j])); off += m[j]
            elif m[j] == dN:
                self.qoff.append((j, None, m[j]))
        self.toff = off if fsp else None
        if fsp:
            off += 1
        self.n = off
        rng = rng or np.random.default_rng(0)
        self.S = rng.normal(size=(nslice, self.n)) + 1j * rng.normal(size=(nslice, self.n))
        self.s0 = rng.normal(size=nslice) + 1j * rng.normal(size=nslice)

    def gpolys(self, z):
        k = self.k; B = z.shape[0]
        c = z[:, 0]; lam = z[:, 1:k]; rr = z[:, k:2 * k - 1]
        gs = [c[:, None].astype(complex)]
        for i in range(k - 1):
            gs.append(np.stack([-lam[:, i] * rr[:, i], lam[:, i]], axis=1))
        return gs

    def F(self, z):
        k = self.k; B = z.shape[0]
        c = z[:, 0]; lam = z[:, 1:k]; rr = z[:, k:2 * k - 1]
        gs = self.gpolys(z)
        A = np.tile(np.array([-1, 0, 1], dtype=complex), (B, 1))
        res = []
        # Korselt
        for j in range(1, k):
            rj = rr[:, j - 1]
            a = rj**2 - 1
            ep = np.ones((B, 1), dtype=complex)
            for i in range(k):
                if i == j:
                    continue
                yi = gs[i][:, 0] + (gs[i][:, 1] * rj if gs[i].shape[1] > 1 else 0)
                ep = pmul(ep, np.stack([np.ones(B), yi], axis=1))
            K = np.zeros(B, dtype=complex)
            for mm in range(k - 1, 0, -1):
                K = K * a + ep[:, mm]
            res.append(K[:, None])
        # Lucas
        for (j, off, mj) in self.qoff:
            N = padd(np.tile(np.array([[2.0]], dtype=complex), (B, 1)), pmul(A, gs[j]))
            pr = np.ones((B, 1), dtype=complex)
            for i in range(k):
                if i != j:
                    pr = pmul(pr, padd(gs[j], -2 * gs[i]))
            gk = np.ones((B, 1), dtype=complex)
            for _ in range(k - 1):
                gk = pmul(gk, gs[j])
            Phi = padd(pr, -gk)
            if off is None:
                lead = N[:, -1:]
                Nm = N / lead
                res.append(prem_monic(Phi, Nm))
            else:
                Q = np.concatenate([z[:, off:off + mj], np.ones((B, 1))], axis=1)
                res.append(prem_monic(N, Q))
                res.append(prem_monic(Phi, Q))
        if self.fsp:
            t = z[:, self.toff]
            for j in range(k):
                gj = gs[j][:, 0] + (gs[j][:, 1] * t if gs[j].shape[1] > 1 else 0)
                res.append((2 + (t**2 - 1) * gj)[:, None])
        if self.nslice:
            res.append(z @ self.S.T - self.s0)
        return np.concatenate(res, axis=1)

    def J(self, z, h=1e-6):
        B, n = z.shape
        out = np.zeros((B, self.M, n), dtype=complex)
        for i in range(n):
            e = np.zeros(n); e[i] = 1
            hh = h * np.maximum(1, np.abs(z[:, i]))[:, None]
            fp = self.F(z + hh * e); fm = self.F(z - hh * e)
            out[:, :, i] = (fp - fm) / (2 * hh)
        return out

def lm_batch(model, z, iters=150, tol=1e-13):
    B, n = z.shape
    with np.errstate(all='ignore'):
        f = model.F(z)
    model.M = f.shape[1]
    nf = np.linalg.norm(f, axis=1)
    active = np.isfinite(nf)
    z[~active] = 0.1; f[~active] = 0; nf[~active] = np.inf
    mu = np.full(B, 1e-3)
    def hess(idx):
        with np.errstate(all='ignore'):
            Jz = model.J(z[idx])
        JH = np.conj(np.transpose(Jz, (0, 2, 1)))
        return JH @ Jz, (JH @ f[idx][:, :, None])[:, :, 0]
    H = np.zeros((B, n, n), dtype=complex); g = np.zeros((B, n), dtype=complex)
    idx0 = np.where(active)[0]
    H[idx0], g[idx0] = hess(idx0)
    eye = np.eye(n)[None]
    for it in range(iters):
        idx = np.where(active & (nf > tol))[0]
        if len(idx) == 0:
            break
        Hm = H[idx] + mu[idx][:, None, None] * eye
        try:
            d = np.linalg.solve(Hm, -g[idx][:, :, None])[:, :, 0]
        except np.linalg.LinAlgError:
            d = np.zeros((len(idx), n), dtype=complex)
            for t, s in enumerate(idx):
                try:
                    d[t] = np.linalg.solve(Hm[t], -g[s])
                except np.linalg.LinAlgError:
                    pass
        zn = z[idx] + d
        with np.errstate(all='ignore'):
            fn = model.F(zn)
            nn = np.linalg.norm(fn, axis=1)
        ok = np.isfinite(nn) & (nn < nf[idx])
        acc = idx[ok]; rej = idx[~ok]
        z[acc] = zn[ok]; f[acc] = fn[ok]; nf[acc] = nn[ok]
        mu[acc] = np.maximum(mu[acc] / 3, 1e-12); mu[rej] *= 4
        if len(acc):
            H[acc], g[acc] = hess(acc)
        active &= (mu < 1e10) & (np.max(np.abs(z), axis=1) < 1e6)
    return z, nf

def gcd_degrees(z, k, tol=1e-7):
    c = z[0]; lam = z[1:k]; rr = z[k:2 * k - 1]
    def g(i, t):
        return c if i == 0 else lam[i - 1] * (t - rr[i - 1])
    out = []; lucas = []
    for j in range(k):
        if j == 0:
            Ncoef = [c, 0, 2 - c]
        else:
            l, r0 = lam[j - 1], rr[j - 1]
            Ncoef = [l, -l * r0, -l, l * r0 + 2]
        roots = np.roots(Ncoef)
        cnt = 0; pts = []
        for t in roots:
            val = np.prod([g(j, t) - 2 * g(i, t) for i in range(k) if i != j]) - g(j, t)**(k - 1)
            scale = abs(g(j, t))**(k - 1) + 1e-300
            if abs(val) / scale < tol:
                cnt += 1; pts.append(t)
        out.append(cnt); lucas.append(pts)
    return out, lucas

def nondeg(z, k, eps=1e-5):
    c = z[0]; lam = z[1:k]; rr = z[k:2 * k - 1]
    if abs(c) < 1e-4 or np.min(np.abs(lam)) < 1e-4 or np.max(np.abs(rr)) > 300:
        return False
    for i in range(k - 1):
        for j in range(i + 1, k - 1):
            if abs(rr[i] - rr[j]) < eps * max(1, abs(rr[i])):
                return False
    return np.max(np.abs(z[:2 * k - 1])) < 1e5

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('k', type=int); ap.add_argument('m'); ap.add_argument('nbatch', type=int)
    ap.add_argument('bs', type=int); ap.add_argument('seed', type=int)
    ap.add_argument('--fsp', action='store_true'); ap.add_argument('--slices', type=int, default=0)
    ap.add_argument('--out', default=None); ap.add_argument('--iters', type=int, default=150)
    a = ap.parse_args()
    k = a.k; m = [int(t) for t in a.m.split(',')]
    rng = np.random.default_rng(a.seed)
    model = Model(k, m, fsp=a.fsp, nslice=a.slices, rng=rng)
    z0 = np.zeros((1, model.n), dtype=complex); model.M = model.F(z0 + 0.3).shape[1]
    print(f"k={k} m={m} fsp={a.fsp} slices={a.slices} unknowns={model.n} equations={model.M} "
          f"expected dim={model.n - model.M}", flush=True)
    found = []
    for b in range(a.nbatch):
        scale = rng.choice([0.5, 1.0, 2.0], size=(a.bs, 1))
        z = (rng.normal(size=(a.bs, model.n)) + 1j * rng.normal(size=(a.bs, model.n))) * scale
        z, nf = lm_batch(model, z, iters=a.iters)
        for s in np.where(nf < 1e-11)[0]:
            zs = z[s]
            if not nondeg(zs, k):
                continue
            Jz = model.J(zs[None])[0]
            sv = np.linalg.svd(Jz, compute_uv=False)
            rank = int(np.sum(sv > 1e-7 * sv[0]))
            gd, pts = gcd_degrees(zs, k)
            found.append((model.n - rank, gd, zs, pts))
        print(f"  batch {b}: converged non-degenerate so far {len(found)}", flush=True)
    cnt = Counter((d, tuple(gd)) for d, gd, zs, pts in found)
    print(f"converged non-degenerate: {len(found)} / {a.nbatch * a.bs}")
    for key, v in sorted(cnt.items()):
        print(f"  local dim {key[0]}  actual gcd degrees {list(key[1])} (sum {sum(key[1])}): {v}")
    if a.out:
        with open(a.out, 'w') as fh:
            for d, gd, zs, pts in found:
                fh.write(repr((d, gd, [complex(t) for t in zs[:2 * k - 1]], [[complex(t) for t in p] for p in pts])) + "\n")

if __name__ == "__main__":
    main()
