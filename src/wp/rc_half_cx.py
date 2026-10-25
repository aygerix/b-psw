"""Complex version of the "-" half-system for cubics in normal form f_j = 2t + 1 + t(t+1)(a_j t + c_j), a_j, c_j in C:
for every j and BOTH roots w of g_j = a_j t^2 + (a_j + c_j) t + (c_j + 2):  F(w) = 1   (2k complex equations, 2k unknowns).
Newton from random complex starts; a DEFLATION factor divides the residual by the distance to the diagonal set
(prod_{i<j} |m_i - m_j|)^(1/npairs), so Newton is repelled from f_1 = ... = f_k. Collect non-diagonal solutions and
report how many are (numerically) real, and whether real ones are valid (a > 0, both quadratics with non-real roots).
Usage: python3 rc_half_cx.py k starts [seed]
"""
import sys, random, cmath, math

def roots2(A, B, C):
    d = cmath.sqrt(B * B - 4 * A * C); return (-B + d) / (2 * A), (-B - d) / (2 * A)

def resid(z, k, deflate=True):
    a = z[:k]; c = z[k:]
    def F(t):
        v = 1 + 0j
        for j in range(k): v *= 2 * t + 1 + t * (t + 1) * (a[j] * t + c[j])
        return v
    R = []
    for j in range(k):
        if abs(a[j]) < 1e-12: return None
        for w in roots2(a[j], a[j] + c[j], c[j] + 2):
            R.append((F(w) - 1) / (1 + abs(w) ** (3 * k - 3)) ** 0.5)
    if deflate:
        prod = 1.0; npairs = 0
        for i in range(k):
            for j in range(i):
                prod *= abs(a[i] - a[j]) + abs(c[i] - c[j]); npairs += 1
        dfl = prod ** (1.0 / npairs)
        R = [r / max(dfl, 1e-300) for r in R]
    return R

def solve(A, b):
    n = len(b); M = [A[i][:] + [b[i]] for i in range(n)]
    for col in range(n):
        p = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[p][col]) < 1e-300: return None
        M[col], M[p] = M[p], M[col]
        for r in range(col + 1, n):
            f = M[r][col] / M[col][col]
            if f:
                for cc in range(col, n + 1): M[r][cc] -= f * M[col][cc]
    x = [0j] * n
    for r in range(n - 1, -1, -1): x[r] = (M[r][n] - sum(M[r][cc] * x[cc] for cc in range(r + 1, n))) / M[r][r]
    return x

def newton(z, k, iters=100):
    for _ in range(iters):
        R = resid(z, k)
        if R is None: return z, False
        nr = max(abs(r) for r in R)
        if nr < 1e-13: return z, True
        J = []
        for col in range(2 * k):
            z2 = z[:]; h = 1e-8 * (1 + abs(z[col])); z2[col] += h
            R2 = resid(z2, k)
            if R2 is None: return z, False
            J.append([(x - y) / h for x, y in zip(R2, R)])
        A = [[J[cc][r] for cc in range(2 * k)] for r in range(2 * k)]
        dz = solve(A, [-r for r in R])
        if dz is None: return z, False
        step = 1.0
        while step > 1e-5:
            z2 = [x + step * d for x, d in zip(z, dz)]; R2 = resid(z2, k)
            if R2 is not None and max(abs(r) for r in R2) < nr: z = z2; break
            step /= 2
        else: return z, False
    R = resid(z, k)
    return z, R is not None and max(abs(r) for r in R) < 1e-10

def main():
    k, starts = int(sys.argv[1]), int(sys.argv[2]); rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    sols = []; conv = 0
    for s in range(starts):
        z = [complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(2 * k)]
        z, ok = newton(z, k)
        if not ok: continue
        raw = resid(z, k, deflate=False)
        if raw is None or max(abs(r) for r in raw) > 1e-9: continue
        conv += 1
        key = tuple(sorted((round(z[j].real, 5), round(z[j].imag, 5), round(z[k + j].real, 5), round(z[k + j].imag, 5)) for j in range(k)))
        if key not in [kk for kk, _ in sols]: sols.append((key, z))
    real = [(kk, z) for kk, z in sols if max(abs(x.imag) for x in z) < 1e-7]
    valid = [(kk, z) for kk, z in real if all(z[j].real > 0 and ((z[k+j].real - z[j].real) ** 2 - 8 * z[j].real) < 0 and (z[k+j].real ** 2 - 8 * z[j].real) < 0 for j in range(k))]
    print(f"k = {k}: {starts} starts, {conv} converged to non-diagonal solutions, {len(sols)} distinct (up to permutation)")
    print(f"  real solutions: {len(real)};  real AND valid (a > 0, special points non-real): {len(valid)}")
    for kk, z in real[:5]: print("   real:", [(round(z[j].real, 5), round(z[k + j].real, 5)) for j in range(k)])
    for kk, z in sols[:3]: print("   sample:", [complex(round(z[j].real, 4), round(z[j].imag, 4)) for j in range(k)], [complex(round(z[k+j].real, 4), round(z[k+j].imag, 4)) for j in range(k)])

if __name__ == "__main__":
    main()
