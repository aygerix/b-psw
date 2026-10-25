"""Complex Newton search for solutions of the k = 5 real-cubic "-" half in Laurent form with FIXED V
(V = y^4 + 15/4 y + 45/16; see laurent_fixed.gp, gen_half_fixedV.py).
Members a_i = 2 + z d_i + z^2 b_i (z = y + 1). S(y) = sum_{r<5} y^(r-1) e_r(a) - e_5(a)(15/4 y + 45/16) has a forced
double zero at y = -1, so S = (y + 1)^2 T(y) with deg T = 9; T == 0 <=> T(y_m) = 0 at 10 distinct points.
Unknowns (b_1..b_5, d_1..d_5) in C^10. Reports converged solutions, flags degenerate ones (some b_i ~ 0, some
alpha_i = a_i(0) ~ 0, or two members equal).
Also checks each solution against the "+" half: by t -> -1 - t, f -> -f the "+" half is the same system at the
transformed members (b_i, -b_i - d_i). A solution of both halves is a (complex) cubic Williams family, k = 5.
Mode "real": real starts (iterates stay real), i.e. a search for REAL solutions of the "-" half.
Usage: python3 hv5_newton.py starts [seed] [real]
"""
import sys, random, cmath

K = 5
PTS = [1.7 * cmath.exp(2j * cmath.pi * (m + 0.25) / 10) for m in range(10)]

def T_at(x, y):
    b, d = x[:K], x[K:]
    z = y + 1
    a = [2 + z * d[i] + z * z * b[i] for i in range(K)]
    # elementary symmetric e_0..e_5 of a
    e = [1 + 0j] + [0j] * K
    for ai in a:
        for r in range(K, 0, -1):
            e[r] += e[r - 1] * ai
    S = sum(y ** (r - 1) * e[r] for r in range(1, K)) - e[K] * (3.75 * y + 2.8125)
    return S / (z * z)

def F(x):
    return [T_at(x, y) for y in PTS]

def solve(A, rhs):
    n = len(rhs); M = [row[:] + [rhs[i]] for i, row in enumerate(A)]
    for c in range(n):
        p = max(range(c, n), key=lambda i: abs(M[i][c]))
        if abs(M[p][c]) < 1e-300: return None
        M[c], M[p] = M[p], M[c]
        for i in range(c + 1, n):
            f = M[i][c] / M[c][c]
            if f != 0:
                for j in range(c, n + 1): M[i][j] -= f * M[c][j]
    out = [0j] * n
    for i in range(n - 1, -1, -1):
        out[i] = (M[i][n] - sum(M[i][j] * out[j] for j in range(i + 1, n))) / M[i][i]
    return out

REAL = False

def newton(x, iters=80):
    for _ in range(iters):
        f = F(x); nf = max(abs(v) for v in f)
        if nf < 1e-11 * (1 + max(abs(v) for v in x)) ** 6: return x, True
        J = [[0j] * 10 for _ in range(10)]
        for c in range(10):
            h = 1e-7 * (1 + abs(x[c]))
            x2 = x[:]; x2[c] += h
            f2 = F(x2)
            for r in range(10): J[r][c] = (f2[r] - f[r]) / h
        dx = solve(J, [-v for v in f])
        if dx is None: return x, False
        # damping
        lam = 1.0
        while lam > 1e-4:
            x2 = [a + lam * b for a, b in zip(x, dx)]
            if REAL: x2 = [complex(v.real, 0) for v in x2]
            if max(abs(v) for v in F(x2)) < nf: break
            lam /= 2
        x = x2
        if max(abs(v) for v in x) > 1e8: return x, False
    f = F(x)
    return x, max(abs(v) for v in f) < 1e-9 * (1 + max(abs(v) for v in x)) ** 6

def plus_resid(x):
    b, d = x[:K], x[K:]
    xt = list(b) + [-b[i] - d[i] for i in range(K)]
    return max(abs(v) for v in F(xt))

def main():
    starts = int(sys.argv[1]); rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    global REAL
    real = len(sys.argv) > 3 and sys.argv[3] == "real"
    REAL = real
    sols = []
    for s in range(starts):
        x = [complex(rng.gauss(0, 2), 0 if real else rng.gauss(0, 2)) for _ in range(10)]
        x, ok = newton(x)
        if not ok: continue
        b, d = x[:K], x[K:]
        alpha = [2 + d[i] + b[i] for i in range(K)]
        key = tuple(sorted((round(b[i].real, 5), round(b[i].imag, 5), round(d[i].real, 5), round(d[i].imag, 5)) for i in range(K)))
        if any(key == s2[0] for s2 in sols): continue
        flags = []
        if min(abs(v) for v in b) < 1e-6: flags.append("b~0")
        if min(abs(v) for v in alpha) < 1e-6: flags.append("alpha~0")
        if any(abs(b[i] - b[j]) + abs(d[i] - d[j]) < 1e-6 for i in range(K) for j in range(i)): flags.append("equal members")
        sols.append((key, flags, b, d))
        isreal = max(abs(v.imag) for v in x) < 1e-8
        print(f"start {s}: solution flags={flags} real={isreal} plus_resid={plus_resid(x):.3e} #{len(sols)}", flush=True)
        print("   b =", [complex(round(v.real, 6), round(v.imag, 6)) for v in b])
        print("   d =", [complex(round(v.real, 6), round(v.imag, 6)) for v in d])
    print(f"{starts} starts: {len(sols)} distinct converged solutions (up to permutation)")

if __name__ == "__main__":
    main()
