"""Numerical search for non-monomial LINEAR Williams families (W'), any odd k. Stdlib only (complex Newton).

Linear members f_j = a_j t + b_j. Write f_j = (t - r_j)/u_j (r_j = root, u_j = 1/a_j), P(t) = prod (t - r_i),
U = prod u_i. Then f_j -+ 1 | F -+ 1 (F = prod f_i), with no coprimality assumed, is equivalent to
    P(r_j + u_j) = U   and   P(r_j - u_j) = -U      for every j        (2k equations).
(Proof: the roots of f_j -+ 1 are r_j +- u_j, and prod_{i != j} f_i(r_j + x) = P(r_j + x) / (x U / u_j) ... = 1.)
Monomial families: all r_j equal, u_j^k = U.
Normalisation: sum r_j = 0 (translation of t), u_1 = 1 (scaling). Unknowns: r_2..r_k, u_2..u_k (2k - 2).
Method: Newton on the square subsystem without the two j = 1 equations, from random complex starts; every
converged point is then checked against the two dropped equations. Reports solutions with r != 0.
Usage: python3 linsearch.py k starts [seed]
"""
import sys, random, cmath

def residuals(k, z):
    r = [0j] * k; u = [1 + 0j] * k
    r[1:] = z[:k - 1]; r[0] = -sum(r[1:])
    u[1:] = z[k - 1:]
    U = 1 + 0j
    for x in u: U *= x
    def P(x):
        v = 1 + 0j
        for ri in r: v *= (x - ri)
        return v
    E = []
    for j in range(k):
        E.append(P(r[j] + u[j]) - U)
        E.append(P(r[j] - u[j]) + U)
    return E, r, u

def solve_lin(A, b):
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
        s = M[i][n] - sum(M[i][jj] * x[jj] for jj in range(i + 1, n))
        x[i] = s / M[i][i]
    return x

def newton(k, z, iters=60):
    m = 2 * k - 2
    for _ in range(iters):
        E, _, _ = residuals(k, z)
        F = E[2:]                       # drop the j = 1 pair
        nf = max(abs(v) for v in F)
        if nf < 1e-13: return z, True
        J = [[0j] * m for _ in range(m)]
        h = 1e-7
        for c in range(m):
            z2 = z[:]; z2[c] += h
            E2, _, _ = residuals(k, z2)
            for i in range(m): J[i][c] = (E2[2 + i] - F[i]) / h
        dz = solve_lin(J, [-v for v in F])
        if dz is None: return z, False
        z = [a + b for a, b in zip(z, dz)]
        if max(abs(v) for v in z) > 1e6: return z, False
    E, _, _ = residuals(k, z)
    return z, max(abs(v) for v in E[2:]) < 1e-10

def main():
    k = int(sys.argv[1]); starts = int(sys.argv[2])
    rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    conv = full = mono = 0; found = []
    for s in range(starts):
        z = [complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(2 * k - 2)]
        z, ok = newton(k, z)
        if not ok: continue
        E, r, u = residuals(k, z)
        if min(abs(x) for x in u) < 1e-4: continue        # degenerate: some u_j = 0, U = 0 (excluded: a_j finite)
        conv += 1
        scale = max(1.0, max(abs(x) for x in u) ** k)
        if max(abs(E[0]), abs(E[1])) / scale > 1e-8: continue
        full += 1
        if max(abs(x) for x in r) < 1e-6:
            mono += 1
        else:
            found.append((r, u))
            if len(found) <= 10:
                print("NON-MONOMIAL candidate: r =", [complex(round(x.real, 6), round(x.imag, 6)) for x in r],
                      " u =", [complex(round(x.real, 6), round(x.imag, 6)) for x in u], flush=True)
    print(f"k = {k}: starts {starts}, non-degenerate square-subsystem roots {conv}, full solutions {full}, monomial {mono}, "
          f"non-monomial {len(found)}")

if __name__ == "__main__":
    main()
