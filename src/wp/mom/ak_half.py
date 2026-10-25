# ak_half.py -- complex Newton for the "-" half of the real-cubic Williams system, general k (all cubic, optional linear member):
#   g_i | (F - 1)/t  for every cubic member i,  f_i = 2t+1 + t(t+1)(a_i t + c_i), g_i = a t^2 + (a+c) t + c + 2.
# 2m complex equations (remainder coefficients) in 2m complex unknowns (a_i, c_i).
# For each converged non-degenerate solution prints the cofactor V = (F-1)/(t prod g_i) (monic, deg k-1) and
# the moments S = sum v, T = sum v u (v = 1/a, u = c/a), and the "+"-half residual.
# usage: python3 ak_half.py m lin starts seed [real]
import sys, numpy as np
m = int(sys.argv[1]); lin = int(sys.argv[2]); N = int(sys.argv[3]); seed = int(sys.argv[4]); REAL = len(sys.argv) > 5
rng = np.random.default_rng(seed); k = m + lin
P = np.polynomial.polynomial   # coefficient arrays, lowest degree first
def members(x):
    a, c = x[:m], x[m:]
    fs = [np.array([1, c[i] + 2, a[i] + c[i], a[i]], dtype=complex) for i in range(m)]
    gs = [np.array([c[i] + 2, a[i] + c[i], a[i]], dtype=complex) for i in range(m)]
    hs = [np.array([2, c[i], a[i]], dtype=complex) for i in range(m)]
    if lin: fs.append(np.array([1, 2], dtype=complex))
    return fs, gs, hs
def prodp(L):
    out = np.array([1], dtype=complex)
    for p in L: out = P.polymul(out, p)
    return out
def remainders(fs, gs, sign):
    # sign = -1: (F-1)/t mod g_i ; sign = +1: (F+1)/(t+1) mod h_i (gs then holds h_i)
    F = prodp(fs)
    if sign < 0:
        Fm = F.copy(); Fm[0] -= 1; R = Fm[1:]                 # (F-1)/t
    else:
        Fp = F.copy(); Fp[0] += 1; R, rr = P.polydiv(Fp, np.array([1, 1], dtype=complex))
    out = []
    for g in gs:
        q, r = P.polydiv(R, g)
        r = np.concatenate([r, np.zeros(2 - len(r))])
        out.extend(list(r))
    return np.array(out)
def resid(x):
    fs, gs, hs = members(x)
    return remainders(fs, gs, -1)
def presid(x):
    fs, gs, hs = members(x)
    return remainders(fs, hs, +1)
def newton(x, iters=100):
    for it in range(iters):
        f = resid(x); nf = np.max(np.abs(f))
        if nf < 1e-12: return x, True
        J = np.zeros((2*m, 2*m), dtype=complex)
        for col in range(2*m):
            h = 1e-7*(1 + abs(x[col])); x2 = x.copy(); x2[col] += h
            J[:, col] = (resid(x2) - f)/h
        try: dx = np.linalg.solve(J, -f)
        except np.linalg.LinAlgError: return x, False
        lam = 1.0
        while lam > 1e-6:
            x2 = x + lam*dx
            if REAL: x2 = x2.real.astype(complex)
            if np.max(np.abs(resid(x2))) < nf: break
            lam /= 2
        x = x2
        if np.max(np.abs(x)) > 1e6: return x, False
    return x, np.max(np.abs(resid(x))) < 1e-10
nsol = 0; Vs = []
for st in range(N):
    if REAL:
        a0 = rng.uniform(0.05, 30, m); c0 = a0/2 + (np.sqrt(8*a0) - a0/2)*rng.uniform(-1, 1, m); x0 = np.concatenate([a0, c0]).astype(complex)
    else:
        x0 = np.concatenate([rng.uniform(0.1, 20, m)*np.exp(1j*rng.uniform(-1, 1, m)), rng.normal(0, 4, m) + 1j*rng.normal(0, 2, m)])
    x, ok = newton(x0)
    if not ok: continue
    a, c = x[:m], x[m:]
    if np.min(np.abs(a)) < 1e-6: continue
    D = np.abs(a[:, None] - a[None, :]) + np.abs(c[:, None] - c[None, :]) + np.eye(m)*1e9
    if D.min() < 1e-5: continue
    fs, gs, hs = members(x)
    F = prodp(fs); Fm = F.copy(); Fm[0] -= 1
    V, rV = P.polydiv(Fm[1:], prodp(gs))
    V = V/V[-1]
    S = np.sum(1/a); T = np.sum(c/a**2)
    nsol += 1
    isreal = np.max(np.abs(x.imag)) < 1e-8
    print("SOL st=%d real=%d  V(low->high)=%s  S=%s T=%s  plusres=%.2e" % (st, isreal, np.array2string(V, precision=6, max_line_width=400), np.round(S, 6), np.round(T, 6), np.max(np.abs(presid(x)))), flush=True)
print("m=%d lin=%d: %d nondegenerate solutions from %d starts" % (m, lin, nsol, N))
