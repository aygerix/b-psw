#!/usr/bin/env python3
"""src/cs/e3b/advchk/sq_enum.py -- adversarial, independent incidence-structure solver with SQUARING-UP (k = 5).

Purpose: the E3b node model (node_model.py) runs Levenberg-Marquardt on OVERDETERMINED systems from random starts; its
recall on known non-empty overdetermined structures is the question.  Here every system is made square (n random
complex combinations of the M >= n equations, Bertini-style randomisation), so that every isolated solution of the true
system is a regular isolated solution of the square system; Newton from random starts then has a positive basin for it.
Every square solution is then checked against ALL M true equations.

Unknowns: c, l_2..l_5, r_2..r_5 (explicit, NOT eliminated), then one position per node.  A = x^2 - 1 (or x^2 with --asq).
Equations: Korselt_j (j = 2..5) in K-form; per node t with member set S:
   first node of member 1: f_1(w) + 1 = 0;   second node of member 1: w' + w = 0;
   first node of cubic j: f_j(th) + 1 = 0;   second: divided difference th^2 + th th' + th'^2 - r_j (th + th') - b = 0;
   third: th + th' + th'' = r_j    (b = 1 for A = x^2 - 1, b = 0 for A = x^2);
   and F(th) + 1 = 0 if |S| < 5.
Optional --pencil "23": impose r_2 = r_3 by substitution (the pencil Korselt component that the E3b nondeg filter drops).
Every converged solution gets a Lucas table recomputed from the polynomials (roots of f_j + 1, |F + 1| there).
Usage: sq_enum.py STRUCT nbatch bs seed [--asq] [--pencil 23] [--out f]"""
import sys, argparse, itertools, numpy as np
from collections import Counter

def parse(s):
    return [tuple(sorted(int(ch) for ch in node if ch.isdigit())) for node in s.split('|')]

class Sys:
    def __init__(self, struct, asq=False, pencil=None, rng=None):
        self.nodes = parse(struct); self.asq = asq; self.b = 0.0 if asq else 1.0
        self.pencil = pencil                      # e.g. (2, 3): r_3 := r_2
        self.N = len(self.nodes)
        self.mn = {j: [t for t, S in enumerate(self.nodes) if j in S] for j in range(1, 6)}
        self.nr = 4 - (1 if pencil else 0)
        self.n = 5 + self.nr + self.N
        M = self.full(np.full((1, self.n), 0.37 + 0.21j)).shape[1]
        self.M = M
        rng = rng or np.random.default_rng(0)
        self.R = None
        if M > self.n:
            self.R = rng.normal(size=(self.n, M)) + 1j * rng.normal(size=(self.n, M))
            self.R[:, :self.n] += 3 * np.eye(self.n)

    def unpack(self, Z):
        c = Z[:, 0]; l = Z[:, 1:5]
        if self.pencil:
            p, q = self.pencil
            rr = []; it = 5
            for j in range(2, 6):
                if j == q:
                    rr.append(None)
                else:
                    rr.append(Z[:, it]); it += 1
            rr[q - 2] = rr[p - 2]
            r = np.stack(rr, axis=1)
        else:
            r = Z[:, 5:9]
        th = Z[:, 5 + self.nr:]
        return c, l, r, th

    def A(self, x):
        return x**2 - self.b

    def g(self, c, l, r, j, x):
        return c + 0 * x if j == 1 else l[:, j - 2] * (x - r[:, j - 2])

    def f(self, c, l, r, j, x):
        return 1 + self.A(x) * self.g(c, l, r, j, x)

    def full(self, Z):
        c, l, r, th = self.unpack(Z)
        res = []
        for j in range(2, 6):
            rj = r[:, j - 2]; a = self.A(rj)
            ys = [self.g(c, l, r, i, rj) for i in range(1, 6) if i != j]
            e = [np.ones_like(c)]
            for y in ys:
                e = [e[0]] + [e[q] + e[q - 1] * y for q in range(1, len(e))] + [e[-1] * y]
            K = np.zeros_like(c)
            for q in range(len(e) - 1, 0, -1):
                K = K * a + e[q]
            res.append(K)
        for t, S in enumerate(self.nodes):
            x = th[:, t]
            for j in S:
                pos = self.mn[j].index(t)
                if pos == 0:
                    res.append(self.f(c, l, r, j, x) + 1)
                elif j == 1:
                    res.append(x + th[:, self.mn[1][0]])
                elif pos == 1:
                    x0 = th[:, self.mn[j][0]]
                    res.append(x0**2 + x0 * x + x**2 - r[:, j - 2] * (x0 + x) - self.b)
                else:
                    x0 = th[:, self.mn[j][0]]; x1 = th[:, self.mn[j][1]]
                    res.append(x0 + x1 + x - r[:, j - 2])
            if len(S) < 5:
                F = np.ones_like(c)
                for i in range(1, 6):
                    F = F * self.f(c, l, r, i, x)
                res.append(F + 1)
        if self.asq:                               # fix the scaling x -> alpha x of A = x^2: r_2 = 1
            res.append(r[:, 0] - 1)
        return np.stack(res, axis=1)

    def sq(self, Z):
        Fz = self.full(Z)
        return Fz if self.R is None else Fz @ self.R.T

    def jac(self, Z, h=1e-7):
        B, n = Z.shape
        cols = []
        for i in range(n):
            e = np.zeros(n); e[i] = 1
            hh = h * np.maximum(1, np.abs(Z[:, i]))[:, None]
            cols.append((self.sq(Z + hh * e) - self.sq(Z - hh * e)) / (2 * hh))
        return np.stack(cols, axis=2)

def newton(S, Z, iters=100):
    with np.errstate(all='ignore'):
        f = S.sq(Z)
    nf = np.linalg.norm(f, axis=1); alive = np.isfinite(nf)
    for it in range(iters):
        idx = np.where(alive & (nf > 1e-13))[0]
        if len(idx) == 0:
            break
        with np.errstate(all='ignore'):
            J = S.jac(Z[idx])
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
            st = np.ones(len(idx)); Zn = Z[idx] + d; fn = S.sq(Zn); nn = np.linalg.norm(fn, axis=1)
            for _ in range(5):
                bad = ~(np.isfinite(nn) & (nn < nf[idx]))
                if not bad.any():
                    break
                st[bad] /= 2; Zn[bad] = Z[idx][bad] + st[bad][:, None] * d[bad]
                fn[bad] = S.sq(Zn[bad]); nn[bad] = np.linalg.norm(fn[bad], axis=1)
        Z[idx] = Zn; f[idx] = fn; nf[idx] = nn
        alive &= np.isfinite(nf) & (np.max(np.abs(Z), axis=1) < 1e4)
        if it == 40:
            alive &= nf < 1e-2
    return Z, nf

def lucas(S, z):
    Z = z[None]; c, l, r, th = S.unpack(Z)
    out = []
    for j in range(1, 6):
        if j == 1:
            coef = [c[0], 0, 2 - S.b * c[0]]
        else:
            lj, rj = l[0, j - 2], r[0, j - 2]
            coef = [lj, -lj * rj, -S.b * lj, S.b * lj * rj + 2]
        row = []
        for x in np.roots(coef):
            xx = np.array([x])
            vals = [S.f(c, l, r, i, xx)[0] for i in range(1, 6)]
            row.append((x, abs(np.prod(vals) + 1) / (np.prod(np.abs(vals)) + 1)))
        out.append(row)
    return out

def degenerate(S, z):
    c, l, r, th = S.unpack(z[None])
    c = c[0]; l = l[0]; r = r[0]
    if abs(c) < 1e-6 or np.min(np.abs(l)) < 1e-6 or np.max(np.abs(z)) > 1e3:
        return 'deg'
    for i in range(4):
        for j in range(i + 1, 4):
            if abs(r[i] - r[j]) < 1e-6 and abs(l[i] - l[j]) < 1e-6:
                return 'deg'
    th = th[0]
    for u in range(len(th)):
        for v in range(u + 1, len(th)):
            if abs(th[u] - th[v]) < 1e-6:
                return 'deg'
    pen = any(abs(r[i] - r[j]) < 1e-6 for i in range(4) for j in range(i + 1, 4))
    return 'pencil' if pen else 'ok'

def famkey(S, z):
    c, l, r, th = S.unpack(z[None]); c = c[0]; l = l[0]; r = r[0]
    def k(sg):
        return tuple([(round((c).real, 5), round((c).imag, 5))] +
                     sorted((round((sg * a).real, 5), round((sg * a).imag, 5), round((sg * b).real, 5), round((sg * b).imag, 5))
                            for a, b in zip(l, r)))
    return min(k(1), k(-1))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('struct'); ap.add_argument('nbatch', type=int); ap.add_argument('bs', type=int); ap.add_argument('seed', type=int)
    ap.add_argument('--asq', action='store_true'); ap.add_argument('--pencil', default=None); ap.add_argument('--out', default=None)
    ap.add_argument('--scale', type=float, default=1.0)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    pen = tuple(int(ch) for ch in a.pencil) if a.pencil else None
    structs = [l.strip() for l in open(a.struct[1:])] if a.struct.startswith('@') else [a.struct]
    for st in structs:
        if not st:
            continue
        S = Sys(st, asq=a.asq, pencil=pen, rng=rng)
        fams = {}; conv = 0; sqsol = 0
        for b in range(a.nbatch):
            Z = (rng.normal(size=(a.bs, S.n)) + 1j * rng.normal(size=(a.bs, S.n))) * rng.choice([0.5, 1, 2], size=(a.bs, 1)) * a.scale
            Z, nf = newton(S, Z)
            for s in np.where(nf < 1e-11)[0]:
                sqsol += 1
                z = Z[s]
                with np.errstate(all='ignore'):
                    tr = np.linalg.norm(S.full(z[None])[0])
                if not (tr < 1e-9):
                    continue
                dg = degenerate(S, z)
                if dg == 'deg':
                    continue
                conv += 1
                key = famkey(S, z)
                if key not in fams:
                    tab = lucas(S, z)
                    gd = [sum(1 for (x, rr) in row if rr < 1e-8) for row in tab]
                    real = np.max(np.abs(z[:5 + S.nr].imag)) < 1e-8
                    fams[key] = (gd, dg, real, z)
        cnt = Counter((tuple(v[0]), v[1]) for v in fams.values())
        top = " ".join(f"{list(k[0])}{'P' if k[1] == 'pencil' else ''}x{v}" for k, v in sorted(cnt.items(), key=lambda kv: -sum(kv[0][0]))[:4])
        best = max((sum(v[0]) for v in fams.values()), default=None)
        print(f"{st:22s} n={S.n} M={S.M} sqsol={sqsol} true={conv} distinct={len(fams)} max_sum={best}  {top}", flush=True)
        if a.out:
            with open(a.out, 'a') as fh:
                for key, (gd, dg, real, z) in fams.items():
                    fh.write(repr((st, gd, dg, real, [complex(w) for w in z])) + "\n")

if __name__ == "__main__":
    main()
