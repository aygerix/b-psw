#!/usr/bin/env python3
"""src/cs/e3b/strat_num.py -- E3b: numerical survey of Lucas gcd-degree strata on shared-A Korselt families.

Family (any odd k): A = x^2 - 1 (WLOG over C when disc A != 0: affine change of x), g_1 = c, g_j = l_j (x - r_j)
(j = 2..k), f_j = 1 + A g_j, F = prod f_j.  Unknowns z = (c, l_2..l_k, r_2..r_k): 2k - 1 complex numbers.
Korselt_j (j >= 2): K_j(A(r_j)) = 0, K_j(a) = [prod_{i != j} (1 + a g_i(r_j)) - 1] / a   (exact: (x - r_j) | (F-1)/A).
Lucas: deg gcd(F + 1, f_j + 1) = deg gcd(N_j, Phi_j), N_j = f_j + 1 = 2 + A g_j, Phi_j = prod_{i != j}(g_j - 2 g_i) - g_j^(k-1)
  [PROVED: g_j^(k-1) (F + 1) == -Phi_j mod N_j and gcd(g_j, N_j) = 1].
Pattern m = (m_1..m_k): gcd degree >= m_j via an auxiliary monic Q_j of degree m_j with N_j == Phi_j == 0 mod Q_j
(m_j = deg N_j: prem(Phi_j, N_j) = 0).  Expected dimension of a stratum: (k - 2)... = dim Korselt - sum m.
Solver: complex Levenberg-Marquardt from random starts; a solution is kept if the residual is < 1e-11 and it is
non-degenerate (c, l_j != 0, r_j distinct, r_j != +-1 not required).  Local dimension = #unknowns - rank J.
Usage: strat_num.py k "m1,m2,...,mk" nstarts seed
"""
import sys, itertools
import numpy as np
import sympy as sp

def build(k, m):
    x = sp.Symbol('x')
    c = sp.Symbol('c')
    lam = sp.symbols('l2:%d' % (k + 1))
    rr = sp.symbols('r2:%d' % (k + 1))
    A = x**2 - 1
    g = [c] + [lam[i] * (x - rr[i]) for i in range(k - 1)]
    unk = [c] + list(lam) + list(rr)
    eqs = []
    a = sp.Symbol('a')
    for j in range(1, k):
        pr = sp.Integer(1)
        for i in range(k):
            if i != j:
                pr *= (1 + a * g[i].subs(x, rr[j - 1]))
        K = sp.expand(sp.cancel((sp.expand(pr) - 1) / a))
        eqs.append(sp.expand(K.subs(a, rr[j - 1]**2 - 1)))
    nkor = len(eqs)
    for j in range(k):
        if m[j] == 0:
            continue
        N = sp.Poly(sp.expand(2 + A * g[j]), x)
        pr = sp.Integer(1)
        for i in range(k):
            if i != j:
                pr *= (g[j] - 2 * g[i])
        Phi = sp.Poly(sp.expand(pr - g[j]**(k - 1)), x)
        dN = N.degree()
        if m[j] > dN:
            raise ValueError("m_j > deg N_j")
        if m[j] == dN:
            R = sp.prem(Phi, N)
            eqs += [sp.expand(cf) for cf in sp.Poly(R, x).all_coeffs()] if not R.is_zero else []
        else:
            qs = sp.symbols('q%d_0:%d' % (j + 1, m[j]))
            unk += list(qs)
            Q = sp.Poly(x**m[j] + sum(qs[t] * x**t for t in range(m[j])), x)
            for P in (N, Phi):
                R = sp.rem(P, Q)            # Q monic: exact polynomial remainder
                if not R.is_zero:
                    eqs += [sp.expand(cf) for cf in R.all_coeffs()]
    F = sp.lambdify(unk, eqs, 'numpy')
    Jm = sp.Matrix(eqs).jacobian(unk)
    J = sp.lambdify(unk, Jm, 'numpy')
    return unk, eqs, F, J, nkor

def lm(F, J, z, iters=200, tol=1e-13):
    mu = 1e-3
    f = np.array(F(*z), dtype=complex)
    nf = np.linalg.norm(f)
    for it in range(iters):
        Jz = np.array(J(*z), dtype=complex)
        H = Jz.conj().T @ Jz
        gvec = Jz.conj().T @ f
        while True:
            try:
                d = np.linalg.solve(H + mu * np.eye(len(z)), -gvec)
            except np.linalg.LinAlgError:
                mu *= 10; continue
            zn = z + d
            fn = np.array(F(*zn), dtype=complex)
            nn = np.linalg.norm(fn)
            if np.isfinite(nn) and nn < nf:
                z, f, nf = zn, fn, nn
                mu = max(mu / 3, 1e-12)
                break
            mu *= 4
            if mu > 1e12:
                return z, nf
        if nf < tol:
            break
    return z, nf

def nondeg(z, k, eps=1e-5):
    c = z[0]; lam = z[1:k]; rr = z[k:2 * k - 1]
    if abs(c) < eps or min(abs(l) for l in lam) < eps:
        return False
    for i in range(k - 1):
        for j in range(i + 1, k - 1):
            if abs(rr[i] - rr[j]) < eps:
                return False
    if np.max(np.abs(z)) > 1e5:
        return False
    return True

def gcd_degrees(z, k, tol=1e-7):
    """Numerical deg gcd(N_j, Phi_j) per member: count roots of N_j where Phi_j vanishes (relative)."""
    c = z[0]; lam = z[1:k]; rr = z[k:2 * k - 1]
    out = []
    for j in range(k):
        def gj(t, j=j):
            return c if j == 0 else lam[j - 1] * (t - rr[j - 1])
        def gi(t, i):
            return c if i == 0 else lam[i - 1] * (t - rr[i - 1])
        if j == 0:
            Ncoef = [c, 0, 2 - c]
        else:
            l, r0 = lam[j - 1], rr[j - 1]
            Ncoef = [l, -l * r0, -l, l * r0 + 2]
        roots = np.roots(Ncoef)
        cnt = 0
        for t in roots:
            val = np.prod([gj(t) - 2 * gi(t, i) for i in range(k) if i != j]) - gj(t)**(k - 1)
            scale = abs(gj(t))**(k - 1) + 1e-300
            if abs(val) / scale < tol:
                cnt += 1
        out.append(cnt)
    return out

def main():
    k = int(sys.argv[1]); m = [int(t) for t in sys.argv[2].split(',')]
    nst = int(sys.argv[3]); seed = int(sys.argv[4])
    rng = np.random.default_rng(seed)
    unk, eqs, F, J, nkor = build(k, m)
    n = len(unk); M = len(eqs)
    print(f"k={k} m={m} unknowns={n} equations={M} expected dim={n - M}  (Korselt dim {2*k-1-nkor}, sum m {sum(m)})", flush=True)
    found = []
    for s in range(nst):
        z0 = (rng.normal(size=n) + 1j * rng.normal(size=n)) * rng.choice([0.5, 1.0, 2.0])
        z, nf = lm(F, J, z0)
        if nf < 1e-11 and nondeg(z, k):
            Jz = np.array(J(*z), dtype=complex)
            sv = np.linalg.svd(Jz, compute_uv=False)
            rank = int(np.sum(sv > 1e-8 * sv[0]))
            gd = gcd_degrees(z, k)
            found.append((n - rank, gd, z))
    print(f"converged non-degenerate: {len(found)} / {nst}")
    from collections import Counter
    cnt = Counter((d, tuple(gd)) for d, gd, z in found)
    for key, v in sorted(cnt.items()):
        print(f"  local dim {key[0]}  actual gcd degrees {list(key[1])} (sum {sum(key[1])}): {v}")
    # print a few representatives
    shown = set()
    for d, gd, z in found:
        key = (d, tuple(gd))
        if key in shown:
            continue
        shown.add(key)
        print("  rep", key, "c,l,r =", np.round(z[:2 * k - 1], 6).tolist())

if __name__ == "__main__":
    main()
