#!/usr/bin/env python3
"""src/cs/e3b/advchk/fsp_tu.py -- adversarial, independent re-derivation of the E3b FSP point counts (k = 5).

Coordinates (T, u_2..u_5) only (no t, no r): scale x = t y.  Then A(ty) = T y^2 - 1 =: a(y), the FSP point is y = 1, and
  g_1 = -2/(T - 1),   g_i(y) = lam_i (y - v_i),  lam_i = -2/((T - 1) u_i),  v_i = 1 - u_i      (i = 2..5)
so every member satisfies f_j(1) = 1 + a(1) g_j(1) = -1.  Korselt_j: K_j(a(v_j)) = 0, K_j(a) = [prod_{i != j}(1 + a g_i(v_j)) - 1]/a.
Extra incidence types:
  '1'  : member 1 Lucas-full, i.e. F(-1) = -1  <=>  prod_i (4 - u_i) = prod_i u_i            (Lemma B; re-derived here
         directly as F(y=-1) + 1 = 0, NOT via that identity)
  'j'  : extra Lucas point s of cubic member j: f_j(s) = -1 with s != 1 (divided difference of N_j at (1, s) = 0) and
         F(s) = -1.
Square systems are solved by damped Newton from random complex starts (batched numpy).  For every converged
non-degenerate solution the FULL Lucas table is recomputed from scratch from the polynomials (roots of f_j + 1 and
|F + 1| there), so Sum m (with distinct-root counting) is independent of the model.
Usage: fsp_tu.py TYPE nbatch bs seed [--out file]     TYPE in {1, 2} (member 2 is WLOG by relabelling)"""
import sys, argparse, numpy as np

def members(Z):
    T = Z[:, 0]; u = Z[:, 1:5]
    g1 = -2 / (T - 1)
    lam = -2 / ((T - 1)[:, None] * u)
    v = 1 - u
    return T, g1, lam, v

def fval(Z, j, y):
    """f_j(t y) in y-coordinates; j = 1..5; y shape (B,)."""
    T, g1, lam, v = members(Z)
    a = T * y**2 - 1
    if j == 1:
        return 1 + a * g1
    return 1 + a * lam[:, j - 2] * (y - v[:, j - 2])

def Fprod(Z, y):
    out = np.ones(Z.shape[0], dtype=complex)
    for j in range(1, 6):
        out = out * fval(Z, j, y)
    return out

def korselt(Z):
    T, g1, lam, v = members(Z)
    res = []
    for j in range(2, 6):
        vj = v[:, j - 2]; a = T * vj**2 - 1
        ys = [g1] + [lam[:, i - 2] * (vj - v[:, i - 2]) for i in range(2, 6) if i != j]
        e = [np.ones_like(T)]                      # elementary symmetric polys of ys
        for y in ys:
            e = [e[0]] + [e[q] + e[q - 1] * y for q in range(1, len(e))] + [e[-1] * y]
        K = np.zeros_like(T)
        for q in range(len(e) - 1, 0, -1):         # K = sum_{q>=1} e_q a^(q-1)
            K = K * a + e[q]
        res.append(K)
    return res

def system(Z, typ):
    res = korselt(Z)
    if typ == '1':
        res.append(Fprod(Z, -np.ones(Z.shape[0], dtype=complex)) + 1)
    else:
        j = int(typ); s = Z[:, 5]
        T, g1, lam, v = members(Z)
        # N_j(y) = 2 + lam (T y^2 - 1)(y - v);  divided difference at (1, s) divided by lam:
        # [ (T s^2 - 1)(s - v) - (T - 1)(1 - v) ] / (s - 1) = T (s^2 + s + 1) - T v (s + 1) - 1
        vj = v[:, j - 2]
        res.append(T * (s**2 + s + 1) - T * vj * (s + 1) - 1)
        res.append(Fprod(Z, s) + 1)
    return np.stack(res, axis=1)

def jac(Z, typ, h=1e-7):
    B, n = Z.shape
    cols = []
    for i in range(n):
        e = np.zeros(n); e[i] = 1
        hh = h * np.maximum(1, np.abs(Z[:, i]))[:, None]
        cols.append((system(Z + hh * e, typ) - system(Z - hh * e, typ)) / (2 * hh))
    return np.stack(cols, axis=2)

def newton(Z, typ, iters=120):
    with np.errstate(all='ignore'):
        f = system(Z, typ)
    nf = np.linalg.norm(f, axis=1)
    alive = np.isfinite(nf)
    for it in range(iters):
        idx = np.where(alive & (nf > 1e-13))[0]
        if len(idx) == 0:
            break
        with np.errstate(all='ignore'):
            J = jac(Z[idx], typ)
            ok = np.all(np.isfinite(J.reshape(len(idx), -1)), axis=1)
            J[~ok] = np.eye(J.shape[1])[None]
            try:
                d = np.linalg.solve(J, -f[idx][:, :, None])[:, :, 0]
            except np.linalg.LinAlgError:
                d = np.zeros((len(idx), Z.shape[1]), dtype=complex)
                for q in range(len(idx)):
                    try:
                        d[q] = np.linalg.solve(J[q], -f[idx[q]])
                    except np.linalg.LinAlgError:
                        pass
            step = np.ones(len(idx))
            # damping: halve until residual decreases (max 6 halvings), else take the small step anyway
            Zn = Z[idx] + d; fn = system(Zn, typ); nn = np.linalg.norm(fn, axis=1)
            for _ in range(6):
                bad = ~(np.isfinite(nn) & (nn < nf[idx]))
                if not bad.any():
                    break
                step[bad] /= 2
                Zn[bad] = Z[idx][bad] + step[bad][:, None] * d[bad]
                fn[bad] = system(Zn[bad], typ); nn[bad] = np.linalg.norm(fn[bad], axis=1)
        Z[idx] = Zn; f[idx] = fn; nf[idx] = nn
        alive &= np.isfinite(nf) & (np.max(np.abs(Z), axis=1) < 1e4)
        if it == 50:
            alive &= nf < 1e-2
    return Z, nf

def lucas_table(Z1):
    """Independent Lucas table in y-coordinates: for each member, roots of f_j + 1 and relative |F + 1| there."""
    Z = Z1[None]
    T, g1, lam, v = members(Z)
    T = T[0]; g1 = g1[0]; lam = lam[0]; v = v[0]
    rows = []
    for j in range(1, 6):
        if j == 1:
            coef = [T * g1, 0, 2 - g1]                      # 2 + g1 (T y^2 - 1)
        else:
            l, vv = lam[j - 2], v[j - 2]
            coef = [l * T, -l * T * vv, -l, l * vv + 2]     # 2 + l (T y^2 - 1)(y - vv)
        rts = np.roots(coef); row = []
        for y in rts:
            val = Fprod(Z, np.array([y]))[0] + 1
            scale = np.prod([abs(fval(Z, i, np.array([y]))[0]) for i in range(1, 6)]) + 1
            row.append((y, abs(val) / scale))
        rows.append(row)
    return rows

def nondeg(z):
    T = z[0]; u = z[1:5]
    if abs(T) < 1e-6 or abs(T - 1) < 1e-6 or np.min(np.abs(u)) < 1e-6 or np.max(np.abs(z[:5])) > 1e3:
        return False
    return min(abs(u[i] - u[j]) for i in range(4) for j in range(i + 1, 4)) > 1e-6

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('typ'); ap.add_argument('nbatch', type=int); ap.add_argument('bs', type=int)
    ap.add_argument('seed', type=int); ap.add_argument('--out', default=None)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    n = 5 if a.typ == '1' else 6
    found = {}
    tot = 0
    for b in range(a.nbatch):
        Z = (rng.normal(size=(a.bs, n)) + 1j * rng.normal(size=(a.bs, n))) * rng.choice([0.3, 1, 3], size=(a.bs, 1))
        Z, nf = newton(Z, a.typ)
        for s in np.where(nf < 1e-11)[0]:
            z = Z[s]
            if not nondeg(z):
                continue
            if a.typ != '1' and abs(z[5] - 1) < 1e-6:
                continue
            tot += 1
            key = (round(z[0].real, 6), round(z[0].imag, 6)) + tuple(sorted((round(w.real, 6), round(w.imag, 6)) for w in z[1:5]))
            if key not in found:
                tab = lucas_table(z)
                gd = [sum(1 for (y, r) in row if r < 1e-8) for row in tab]
                near = sorted(r for row in tab for (y, r) in row if r >= 1e-8)[:1]
                J = jac(z[None], a.typ)[0]
                sv = np.linalg.svd(J, compute_uv=False)
                found[key] = (gd, z, near[0] if near else None, sv[-1] / sv[0])
        print(f"batch {b}: converged non-degenerate {tot}, distinct {len(found)}", flush=True)
    from collections import Counter
    cnt = Counter(tuple(v[0]) for v in found.values())
    print("distinct families by Lucas gcd degrees (y-chart; relabelling of 2..5 identified):")
    for gd, c in sorted(cnt.items(), key=lambda kv: -sum(kv[0])):
        print(f"  {list(gd)} sum {sum(gd)}: {c}")
    nreal = sum(1 for v in found.values() if np.max(np.abs(v[1][:5].imag)) < 1e-8)
    print("real families:", nreal)
    print("max Sum m:", max(sum(v[0]) for v in found.values()) if found else None)
    print("smallest non-Lucas relative residual over all points:", min((v[2] for v in found.values() if v[2] is not None), default=None))
    print("smallest relative singular value of the Jacobian:", min((v[3] for v in found.values()), default=None))
    if a.out:
        with open(a.out, 'w') as fh:
            for gd, z, near, sv in found.values():
                fh.write(repr((gd, [complex(w) for w in z], near, sv)) + "\n")

if __name__ == "__main__":
    main()
