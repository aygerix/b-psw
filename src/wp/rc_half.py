"""Two-stage numerical study of REAL CUBIC Williams families (normal form f = 2t + 1 + t(t+1)(a t + c)).
Stage 1: solve the SQUARE "-" half-system  F(w_j) = 1  (j = 1..k; w_j the upper "-" point of member j) by Newton from
random starts, parametrising member j by w_j = x_j + i exp(y_j) (2k real unknowns, 2k real equations).
Note: the diagonal (all members equal) solves the "-" half; solutions are kept only if members are pairwise distinct
(min separation > SEP) and valid (a_j > 0, "+" points non-real).
Stage 2: for each distinct valid solution, the "+" residual R+ = max_j |F(v_j) + 1| (zero for a Williams family).
Reports: number of distinct real half-solutions (up to permutation, rounded) and the smallest R+.
Usage: python3 rc_half.py k starts [seed]
"""
import sys, random, math, cmath

SEP = 1e-3

def members(z, k):
    out = []
    for j in range(k):
        w = complex(z[2 * j], math.exp(max(-30.0, min(8.0, z[2 * j + 1]))))
        lam = -2 / (w + 1); a = lam.imag / w.imag; c = lam.real - a * w.real
        out.append((a, c, w))
    return out

def F(m, t):
    v = 1 + 0j
    for a, c, _ in m: v *= 2 * t + 1 + t * (t + 1) * (a * t + c)
    return v

def half(z, k):
    m = members(z, k); R = []
    for a, c, w in m:
        e = F(m, w) - 1; s = 1.0 / (1 + abs(w) ** (3 * k - 3)) ** 0.5
        R += [e.real * s, e.imag * s]
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
    x = [0.0] * n
    for r in range(n - 1, -1, -1): x[r] = (M[r][n] - sum(M[r][cc] * x[cc] for cc in range(r + 1, n))) / M[r][r]
    return x

def newton(z, k, iters=80):
    for _ in range(iters):
        R = half(z, k); nr = max(abs(r) for r in R)
        if nr < 1e-14: return z, True
        J = []
        for col in range(2 * k):
            z2 = z[:]; h = 1e-8 * (1 + abs(z[col])); z2[col] += h
            R2 = half(z2, k); J.append([(x - y) / h for x, y in zip(R2, R)])
        A = [[J[c][r] for c in range(2 * k)] for r in range(2 * k)]
        dz = solve(A, [-r for r in R])
        if dz is None: return z, False
        step = 1.0
        while step > 1e-4:
            z2 = [x + step * d for x, d in zip(z, dz)]
            if max(abs(r) for r in half(z2, k)) < nr: z = z2; break
            step /= 2
        else: return z, False
    return z, max(abs(r) for r in half(z, k)) < 1e-11

def main():
    k, starts = int(sys.argv[1]), int(sys.argv[2]); rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    sols = {}; conv = diag = invalid = 0
    for s in range(starts):
        z = []
        for j in range(k): z += [rng.uniform(-3, 2), rng.uniform(-2, 2)]
        z, ok = newton(z, k)
        if not ok: continue
        conv += 1
        m = members(z, k)
        if min(abs(m[i][0] - m[j][0]) + abs(m[i][1] - m[j][1]) for i in range(k) for j in range(i)) < SEP: diag += 1; continue
        if any(a <= 0 or c * c - 8 * a >= 0 for a, c, _ in m): invalid += 1; continue
        Rp = 0.0
        for a, c, w in m:
            v = complex(-c / (2 * a), math.sqrt(8 * a - c * c) / (2 * a))
            Rp = max(Rp, abs(F(m, v) + 1) / (1 + abs(v) ** (3 * k - 3)) ** 0.5)
        key = tuple(sorted((round(a, 6), round(c, 6)) for a, c, _ in m))
        if key not in sols: sols[key] = Rp
    print(f"k = {k}: {starts} starts, {conv} converged, {diag} near-diagonal, {invalid} invalid, {len(sols)} distinct valid real half-solutions")
    for key, Rp in sorted(sols.items(), key=lambda kv: kv[1])[:8]:
        print(f"  R+ = {Rp:.3e}  members (a, c) = {list(key)}")
    print("Williams families found (R+ < 1e-10):", sum(1 for v in sols.values() if v < 1e-10))

if __name__ == "__main__":
    main()
