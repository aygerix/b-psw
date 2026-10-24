"""Search for NON-monomial linear Williams families, any odd k (complex Levenberg-Marquardt, stdlib only).

Same equations as linsearch.py: with P(t) = prod (t - r_i) and U = prod u_i,
    P(r_j + u_j) = U,  P(r_j - u_j) = -U   (j = 1..k), u_j != 0.
Monomial families have all r_j equal. Chart: sum r_j = 0 and prod u_j = 1 (translation and scaling of t). Unknowns r_2..r_k, u_2..u_k (2k - 2);
equations 2k. End points with zero residual and r = 0 are monomial and are counted separately.
Every start runs damped Gauss-Newton (LM) on the overdetermined system. A zero residual is a non-monomial
family, i.e. a counterexample to W' in the linear case. Prints the smallest residuals found.
Usage: python3 linsearch2.py k starts [seed]
"""
import sys, random

def unpack(k, z):
        # chart: sum r_j = 0 (translation), prod u_j = 1 (scaling; keeps away from the degenerate locus U = 0)
    rest = z[:k - 1]
    r = [-sum(rest)] + rest
    uu = z[k - 1:]
    pr = 1 + 0j
    for x in uu: pr *= x
    u = uu + [1 / pr if abs(pr) > 1e-300 else 1e300 + 0j]
    return r, u

def resid(k, z):
    r, u = unpack(k, z)
    U = 1 + 0j
    for x in u: U *= x
    def P(x):
        v = 1 + 0j
        for ri in r: v *= (x - ri)
        return v
    sc = 1.0
    E = []
    for j in range(k):
        E.append((P(r[j] + u[j]) - U) * sc)
        E.append((P(r[j] - u[j]) + U) * sc)
    # keep u away from 0: penalise tiny u (the degenerate locus U = 0)
    return E, r, u, U

def solve(A, b):
    n = len(b); M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        p = max(range(c, n), key=lambda i: abs(M[i][c]))
        if abs(M[p][c]) < 1e-300: return None
        M[c], M[p] = M[p], M[c]
        for i in range(c + 1, n):
            f = M[i][c] / M[c][c]
            if f != 0:
                for jj in range(c, n + 1): M[i][jj] -= f * M[c][jj]
    x = [0j] * n
    for i in range(n - 1, -1, -1):
        x[i] = (M[i][n] - sum(M[i][jj] * x[jj] for jj in range(i + 1, n))) / M[i][i]
    return x

def lm(k, z, iters=200):
    m = 2 * k - 2; lam = 1e-3
    E, _, _, _ = resid(k, z); f = sum(abs(e) ** 2 for e in E)
    for _ in range(iters):
        h = 1e-7
        J = []
        for c in range(m):
            z2 = z[:]; z2[c] += h
            E2, _, _, _ = resid(k, z2)
            J.append([(a - b) / h for a, b in zip(E2, E)])      # column c
        # normal equations (J^H J + lam I) dz = -J^H E
        A = [[sum(J[a][i].conjugate() * J[b][i] for i in range(len(E))) for b in range(m)] for a in range(m)]
        g = [-sum(J[a][i].conjugate() * E[i] for i in range(len(E))) for a in range(m)]
        while True:
            A2 = [row[:] for row in A]
            for a in range(m): A2[a][a] += lam * (1 + abs(A[a][a]))
            dz = solve(A2, g)
            if dz is None: return z, f
            z2 = [a + b for a, b in zip(z, dz)]
            E2, r2, u2, U2 = resid(k, z2); f2 = sum(abs(e) ** 2 for e in E2)
            if f2 < f:
                z, E, f = z2, E2, f2; lam = max(lam / 3, 1e-12); break
            lam *= 4
            if lam > 1e8: return z, f
        if f < 1e-28: break
    return z, f

def main():
    k = int(sys.argv[1]); starts = int(sys.argv[2])
    rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    best = []; mono_hits = [0]
    for s in range(starts):
        z = [complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(2 * k - 2)]
        z, f = lm(k, z)
        r, u = unpack(k, z)
        if min(abs(x) for x in u) < 1e-3: continue                # degenerate (U = 0) locus
        if f < 1e-20 and max(abs(x) for x in r) < 1e-6:
            mono_hits[0] += 1; continue                           # monomial solution
        best.append((f, r, u))
    best.sort(key=lambda t: t[0])
    print(f"k = {k}: {starts} starts, {len(best)} non-degenerate end points; smallest residuals^2:")
    for f, r, u in best[:5]:
        print(f"  {f:.3e}  r = {[complex(round(x.real, 4), round(x.imag, 4)) for x in r]}")
        print(f"           u = {[complex(round(x.real, 4), round(x.imag, 4)) for x in u]}")
    print(f"monomial solutions reached: {mono_hits[0]}")
    zeros = sum(1 for f, _, _ in best if f < 1e-20)
    print(f"end points with residual^2 < 1e-20 (non-monomial solutions): {zeros}")

if __name__ == "__main__":
    main()
