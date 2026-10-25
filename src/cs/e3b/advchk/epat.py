#!/usr/bin/env python3
"""src/cs/e3b/advchk/epat.py -- adversarial gcd-pattern solver WITHOUT auxiliary unknowns (k = 5, A = x^2 - 1).

For member j let C_j be the companion matrix of the monic version of N_j = f_j + 1 and M_j = Phi_j(C_j), whose eigenvalues
are Phi_j(theta) over the roots theta of N_j (Lemma A: theta is a Lucas point iff Phi_j(theta) = 0).  With e_1, e_2, e_3 the
elementary symmetric functions of those eigenvalues (trace / principal minors / det, all polynomial in the parameters):
   m_j >= 1  <=>  e_top = 0;   m_j >= 2 (cubic) <=> e_3 = e_2 = 0;   m_j = 3 <=> e_1 = e_2 = e_3 = 0;  member 1: m_1 = 2 <=> e_1 = e_2 = 0.
So a pattern with Sum m = 5 on the 5-dim Korselt variety is a SQUARE system in the 9 parameters (c, l_2..l_5, r_2..r_5).
Loci (--locus): 'main' (9 params, 4 Korselt eqs); 'pen23' (r_3 = r_2: 8 params, K_2 = K_3 so 3 Korselt eqs);
'symdp' (reflection-symmetric double pencil: params c, l_2, l_3, r; l_4 = -l_2, l_5 = -l_3, r_2 = r_3 = r, r_4 = r_5 = -r;
1 Korselt eq; patterns on members 1, 2, 3 only, mirrors automatic); 'sym' (reflection-symmetric, distinct r: params c, l_2, r_2,
l_4, r_4; f_3 = f_2(-x), f_5 = f_4(-x); 2 Korselt eqs; patterns on members 1, 2, 4).
Every converged solution: full Lucas table recomputed from the polynomials; Sum m recorded; families deduplicated.
Usage: epat.py LOCUS "m1,m2,..." nbatch bs seed [--out f]   (pattern lists the independent members of the locus)"""
import sys, argparse, numpy as np
from collections import Counter

LOCI = {
    'main': dict(npar=9, kor=[2, 3, 4, 5], mem=[1, 2, 3, 4, 5]),
    'pen23': dict(npar=8, kor=[2, 4, 5], mem=[1, 2, 3, 4, 5]),
    'symdp': dict(npar=4, kor=[2], mem=[1, 2, 3]),
    'sym': dict(npar=5, kor=[2, 4], mem=[1, 2, 4]),
}

def expand(locus, P):
    """parameters -> (c, l[B,4], r[B,4])"""
    if locus == 'main':
        return P[:, 0], P[:, 1:5], P[:, 5:9]
    if locus == 'pen23':
        r = np.stack([P[:, 5], P[:, 5], P[:, 6], P[:, 7]], axis=1)
        return P[:, 0], P[:, 1:5], r
    if locus == 'symdp':
        c, l2, l3, r = P[:, 0], P[:, 1], P[:, 2], P[:, 3]
        return c, np.stack([l2, l3, -l2, -l3], axis=1), np.stack([r, r, -r, -r], axis=1)
    if locus == 'sym':
        c, l2, r2, l4, r4 = P[:, 0], P[:, 1], P[:, 2], P[:, 3], P[:, 4]
        return c, np.stack([l2, -l2, l4, -l4], axis=1), np.stack([r2, -r2, r4, -r4], axis=1)

def pmul(a, b):
    out = np.zeros((a.shape[0], a.shape[1] + b.shape[1] - 1), dtype=complex)
    for i in range(a.shape[1]):
        out[:, i:i + b.shape[1]] += a[:, i:i + 1] * b
    return out

def gpoly(c, l, r, j):
    B = c.shape[0]
    if j == 1:
        return c[:, None].astype(complex)
    return np.stack([-l[:, j - 2] * r[:, j - 2], l[:, j - 2]], axis=1).astype(complex)

def padd(a, b, sb=1.0):
    n = max(a.shape[1], b.shape[1]); out = np.zeros((a.shape[0], n), dtype=complex)
    out[:, :a.shape[1]] += a; out[:, :b.shape[1]] += sb * b
    return out

def phi(c, l, r, j):
    gj = gpoly(c, l, r, j)
    pr = np.ones((c.shape[0], 1), dtype=complex)
    for i in range(1, 6):
        if i != j:
            pr = pmul(pr, padd(gj, gpoly(c, l, r, i), -2.0))
    gk = np.ones((c.shape[0], 1), dtype=complex)
    for _ in range(4):
        gk = pmul(gk, gj)
    return padd(pr, gk, -1.0)            # coefficients, lowest degree first

def companion(c, l, r, j):
    B = c.shape[0]
    if j == 1:
        a0 = (2 - c) / c                  # x^2 + a0
        C = np.zeros((B, 2, 2), dtype=complex); C[:, 0, 1] = -a0; C[:, 1, 0] = 1
        return C
    lj, rj = l[:, j - 2], r[:, j - 2]
    # x^3 - r x^2 - x + (r + 2/l)
    a = [rj + 2 / lj, -np.ones(B), -rj]   # a0, a1, a2
    C = np.zeros((B, 3, 3), dtype=complex)
    C[:, 1, 0] = 1; C[:, 2, 1] = 1
    C[:, 0, 2] = -a[0]; C[:, 1, 2] = -a[1]; C[:, 2, 2] = -a[2]
    return C

def esym(c, l, r, j):
    C = companion(c, l, r, j); P = phi(c, l, r, j)
    n = C.shape[1]; I = np.broadcast_to(np.eye(n), C.shape)
    M = np.zeros_like(C)
    for t in range(P.shape[1] - 1, -1, -1):
        M = M @ C + P[:, t][:, None, None] * I
    tr = np.trace(M, axis1=1, axis2=2)
    tr2 = np.trace(M @ M, axis1=1, axis2=2)
    e2 = (tr**2 - tr2) / 2
    det = np.linalg.det(M)
    return [tr, e2, det] if n == 3 else [tr, det]

def korselt(c, l, r, j):
    rj = r[:, j - 2]; a = rj**2 - 1
    ys = [c] + [l[:, i - 2] * (rj - r[:, i - 2]) for i in range(2, 6) if i != j]
    e = [np.ones_like(c)]
    for y in ys:
        e = [e[0]] + [e[q] + e[q - 1] * y for q in range(1, len(e))] + [e[-1] * y]
    K = np.zeros_like(c)
    for q in range(len(e) - 1, 0, -1):
        K = K * a + e[q]
    return K

class Pat:
    def __init__(self, locus, m):
        self.locus = locus; self.L = LOCI[locus]; self.m = dict(zip(self.L['mem'], m)); self.n = self.L['npar']
    def F(self, P):
        c, l, r = expand(self.locus, P)
        res = [korselt(c, l, r, j) for j in self.L['kor']]
        for j, mj in self.m.items():
            if mj == 0:
                continue
            e = esym(c, l, r, j)
            if j == 1:
                if self.locus in ('symdp', 'sym'):      # Phi_1 even here: both roots of N_1 carry the same value
                    res += [e[0]] if mj == 2 else []
                else:
                    res += [e[1]] if mj == 1 else [e[0], e[1]]
            else:
                res += {1: [e[2]], 2: [e[2], e[1]], 3: [e[2], e[1], e[0]]}[mj]
        return np.stack(res, axis=1)
    def J(self, P, h=1e-7):
        cols = []
        for i in range(self.n):
            e = np.zeros(self.n); e[i] = 1
            hh = h * np.maximum(1, np.abs(P[:, i]))[:, None]
            cols.append((self.F(P + hh * e) - self.F(P - hh * e)) / (2 * hh))
        return np.stack(cols, axis=2)

def newton(S, P, iters=100):
    with np.errstate(all='ignore'):
        f = S.F(P)
    nf = np.linalg.norm(f, axis=1); alive = np.isfinite(nf)
    for it in range(iters):
        idx = np.where(alive & (nf > 1e-14 * (1 + 0 * nf)))[0]
        if len(idx) == 0:
            break
        with np.errstate(all='ignore'):
            J = S.J(P[idx])
            okJ = np.all(np.isfinite(J.reshape(len(idx), -1)), axis=1)
            J[~okJ] = np.eye(S.n)[None]
            d = np.zeros((len(idx), S.n), dtype=complex)
            try:
                d = np.linalg.solve(J, -f[idx][:, :, None])[:, :, 0]
            except np.linalg.LinAlgError:
                for q in range(len(idx)):
                    try:
                        d[q] = np.linalg.solve(J[q], -f[idx[q]])
                    except np.linalg.LinAlgError:
                        pass
            st = np.ones(len(idx)); Pn = P[idx] + d; fn = S.F(Pn); nn = np.linalg.norm(fn, axis=1)
            for _ in range(5):
                bad = ~(np.isfinite(nn) & (nn < nf[idx]))
                if not bad.any():
                    break
                st[bad] /= 2; Pn[bad] = P[idx][bad] + st[bad][:, None] * d[bad]
                fn[bad] = S.F(Pn[bad]); nn[bad] = np.linalg.norm(fn[bad], axis=1)
        P[idx] = Pn; f[idx] = fn; nf[idx] = nn
        alive &= np.isfinite(nf) & (np.max(np.abs(P), axis=1) < 1e4)
        if it == 40:
            alive &= nf < 1e-2
    return P, nf

def lucas_sum(c, l, r):
    A = lambda x: x**2 - 1
    def f(j, x):
        return 1 + A(x) * (c if j == 1 else l[j - 2] * (x - r[j - 2]))
    gd = []
    for j in range(1, 6):
        coef = [c, 0, 2 - c] if j == 1 else [l[j - 2], -l[j - 2] * r[j - 2], -l[j - 2], l[j - 2] * r[j - 2] + 2]
        cnt = 0
        for x in np.roots(coef):
            vals = [f(i, x) for i in range(1, 6)]
            if abs(np.prod(vals) + 1) / (np.prod(np.abs(vals)) + 1) < 1e-8:
                cnt += 1
        gd.append(cnt)
    return gd

def nondeg(c, l, r):
    if abs(c) < 1e-5 or np.min(np.abs(l)) < 1e-5 or np.max(np.abs(r)) > 1e3 or np.max(np.abs(l)) > 1e4:
        return False
    for i in range(4):
        for j in range(i + 1, 4):
            if abs(r[i] - r[j]) < 1e-5 and abs(l[i] - l[j]) < 1e-5 * max(1, abs(l[i])):
                return False
    return True

def key(c, l, r):
    def k(s):
        return tuple([(round(c.real, 5), round(c.imag, 5))] + sorted((round((s * a).real, 5), round((s * a).imag, 5),
                       round((s * b).real, 5), round((s * b).imag, 5)) for a, b in zip(l, r)))
    return min(k(1), k(-1))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('locus'); ap.add_argument('m'); ap.add_argument('nbatch', type=int); ap.add_argument('bs', type=int)
    ap.add_argument('seed', type=int); ap.add_argument('--out', default=None)
    a = ap.parse_args()
    m = [int(t) for t in a.m.split(',')]
    S = Pat(a.locus, m)
    rng = np.random.default_rng(a.seed)
    M = S.F(np.full((1, S.n), 0.31 + 0.17j)).shape[1]
    fams = {}; conv = 0
    for b in range(a.nbatch):
        P = (rng.normal(size=(a.bs, S.n)) + 1j * rng.normal(size=(a.bs, S.n))) * rng.choice([0.5, 1, 2], size=(a.bs, 1))
        P, nf = newton(S, P)
        for s in np.where(nf < 1e-10)[0]:
            c, l, r = expand(a.locus, P[s][None]); c = c[0]; l = l[0]; r = r[0]
            if not nondeg(c, l, r):
                continue
            conv += 1
            kk = key(c, l, r)
            if kk not in fams:
                gd = lucas_sum(c, l, r)
                fams[kk] = (gd, c, l, r)
    cnt = Counter(tuple(v[0]) for v in fams.values())
    best = max((sum(v[0]) for v in fams.values()), default=None)
    pen = sum(1 for v in fams.values() if min(abs(v[3][i] - v[3][j]) for i in range(4) for j in range(i + 1, 4)) < 1e-5)
    print(f"{a.locus} m={m} n={S.n} M={M} conv={conv}/{a.nbatch * a.bs} distinct={len(fams)} (pencil-type {pen}) max_sum={best}  " +
          " ".join(f"{list(k)}x{v}" for k, v in sorted(cnt.items(), key=lambda kv: -sum(kv[0]))[:5]), flush=True)
    if a.out:
        with open(a.out, 'a') as fh:
            for gd, c, l, r in fams.values():
                fh.write(repr((a.locus, m, gd, complex(c), [complex(v) for v in l], [complex(v) for v in r])) + "\n")

if __name__ == "__main__":
    main()
