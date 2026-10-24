"""Real CUBIC Williams families in the Lemma 24 normal form, numerical search (real LM, stdlib only).
Members f_j = 2t + 1 + t(t+1)(a_j t + c_j), a_j > 0 (shared real special points t = 0 (f = 1), t = -1 (f = -1)).
Parametrisation (v2): each member by its upper "-" point w_j = x_j + i exp(y_j) (non-real by construction);
(w+1)(a w + c) = -2 gives a = Im(lam)/Im(w), c = Re(lam) - a Re(w), lam = -2/(w+1).
"+" points v_j: roots of a v^2 + c v + 2 (upper one). Residuals F(w_j) - 1, F(v_j) + 1 (real, imag): 4k equations,
2k unknowns. Penalties: a_j > 0; the "+" discriminant c^2 - 8a must stay clearly negative (relative margin 1e-2), and
Im w_j >= 1e-2 (1 + |w_j|), so that LM cannot converge to the real-double-root boundary where the value condition
alone is not the Williams condition (v1 artifact: k = 3 "zeros" with disc g = 0); distinct members.
Usage: python3 rc_lm.py k starts [seed]
"""
import sys, random, cmath, math

def members(z, k):
    out = []
    for j in range(k):
        w = complex(z[2 * j], math.exp(max(-30.0, min(10.0, z[2 * j + 1]))))
        lam = -2 / (w + 1)
        a = lam.imag / w.imag; c = lam.real - a * w.real
        out.append((a, c, w))
    return out

def unpack(z, k):
    m = members(z, k)
    return [x[0] for x in m], [x[1] for x in m]

def resid(z, k):
    m = members(z, k)
    a = [x[0] for x in m]; c = [x[1] for x in m]
    def F(t):
        v = 1 + 0j
        for j in range(k): v *= 2 * t + 1 + t * (t + 1) * (a[j] * t + c[j])
        return v
    R = []
    for j in range(k):
        aj, cj, w = m[j]
        pen = 0.0
        if aj <= 1e-6: pen += (1e-6 - aj) * 1e3 + 1.0
        d2 = cj * cj - 8 * aj
        pen += max(0.0, d2 / (cj * cj + 8 * abs(aj) + 1e-300) + 1e-2) * 10
        pen += max(0.0, 1e-2 * (1 + abs(w)) - w.imag) * 10
        if aj <= 1e-6 or d2 >= 0:
            R += [pen + 1.0, pen + 1.0, pen + 1.0, pen + 1.0, pen]; continue
        v = complex(-cj / (2 * aj), math.sqrt(-d2) / (2 * aj))
        e1 = F(w) - 1; e2 = F(v) + 1
        s = 1.0 / (1 + abs(w) ** (3 * k - 3) + abs(v) ** (3 * k - 3)) ** 0.5
        R += [e1.real * s, e1.imag * s, e2.real * s, e2.imag * s, pen]
    for i in range(k):
        for j in range(i):
            sep = abs(a[i] - a[j]) + abs(c[i] - c[j])
            R.append(max(0.0, 0.05 - sep) * 10)
    return R

def solve(M, b):
    n = len(b); A = [M[i][:] + [b[i]] for i in range(n)]
    for col in range(n):
        p = max(range(col, n), key=lambda r: abs(A[r][col]))
        if abs(A[p][col]) < 1e-300: return None
        A[col], A[p] = A[p], A[col]
        for r in range(col + 1, n):
            f = A[r][col] / A[col][col]
            if f:
                for cc in range(col, n + 1): A[r][cc] -= f * A[col][cc]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        x[r] = (A[r][n] - sum(A[r][cc] * x[cc] for cc in range(r + 1, n))) / A[r][r]
    return x

def lm(z, k, iters=300):
    R = resid(z, k); f = sum(r * r for r in R); lam = 1e-3; m = len(z)
    for _ in range(iters):
        J = []
        for col in range(m):
            z2 = z[:]; h = 1e-7 * (1 + abs(z[col])); z2[col] += h
            R2 = resid(z2, k); J.append([(x - y) / h for x, y in zip(R2, R)])
        A = [[sum(J[p][i] * J[q][i] for i in range(len(R))) for q in range(m)] for p in range(m)]
        g = [-sum(J[p][i] * R[i] for i in range(len(R))) for p in range(m)]
        while True:
            A2 = [row[:] for row in A]
            for p in range(m): A2[p][p] += lam * (1 + A[p][p])
            dz = solve(A2, g)
            if dz is None: return z, f
            z2 = [x + y for x, y in zip(z, dz)]
            R2 = resid(z2, k); f2 = sum(r * r for r in R2)
            if f2 < f: z, R, f = z2, R2, f2; lam = max(lam / 3, 1e-12); break
            lam *= 4
            if lam > 1e10: return z, f
        if f < 1e-26: break
    return z, f

def main():
    k, starts = int(sys.argv[1]), int(sys.argv[2]); rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    res = []
    for s in range(starts):
        z = []
        for j in range(k):
            z += [rng.uniform(-2, 1), rng.uniform(-1.5, 1.5)]
        z, f = lm(z, k)
        res.append((f, z))
    res.sort(key=lambda t: t[0])
    print(f"k = {k}: {starts} starts; smallest residual^2:")
    for f, z in res[:6]:
        a, c = unpack(z, k)
        print(f"  {f:.3e}  a = {[round(x, 4) for x in a]}  c = {[round(x, 4) for x in c]}")
    print(f"zeros (residual^2 < 1e-20): {sum(1 for f, _ in res if f < 1e-20)}")

if __name__ == "__main__":
    main()
