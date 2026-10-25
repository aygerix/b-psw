# ak_half2.py -- "-" half in Laurent form: [s^n] prod_i (1 + x_i(s)) = 0 for n = k..3k-1,
#   x_i = v_i s^3 / (1 + (1+u_i) s + (u_i + 2 v_i) s^2)   (v = 1/a, u = c/a); linear member x = s/2.
# Equivalent to g_i | (F-1)/t for all i (2k equations, 2k unknowns (u_i, v_i)).
# Prints, per non-degenerate solution: V coefficients (= G_0..G_{k-1} as t^{k-1-n}), S, T, and the "+"-half residual
# (the "+" half = the same equations at u -> 1 - u).
# usage: python3 ak_half2.py m lin starts seed [real]
import sys, numpy as np
m = int(sys.argv[1]); lin = int(sys.argv[2]); N = int(sys.argv[3]); seed = int(sys.argv[4]); REAL = len(sys.argv) > 5
rng = np.random.default_rng(seed); k = m + lin; L = 3*k
def xser(u, v):
    c = np.zeros(L, dtype=complex); c[0] = 1; c[1] = -(1 + u)
    for j in range(2, L): c[j] = -(1 + u)*c[j-1] - (u + 2*v)*c[j-2]
    x = np.zeros(L, dtype=complex); x[3:] = v*c[:L-3]; return x
def Gser(U, W):
    G = np.zeros(L, dtype=complex); G[0] = 1
    for i in range(m):
        x = xser(U[i], W[i]); x[0] += 1
        G = np.convolve(G, x)[:L]
    if lin:
        G = np.convolve(G, np.array([1, 0.5]))[:L]
    return G
def resid(z):
    U, W = z[:m], z[m:]
    return Gser(U, W)[k:]
def dfac(z):
    U, W = z[:m], z[m:]
    f = np.prod(1 + 0.05/W)
    for i in range(m):
        for l in range(i+1, m):
            f *= 1 + 0.05/((U[i]-U[l])**2 + (W[i]-W[l])**2)
    return f
def dresid(z):
    return resid(z)*dfac(z)
def presid(z):
    U, W = z[:m], z[m:]
    return Gser(1 - U, W)[k:]
def newton(z, iters=150):
    for it in range(iters):
        f = dresid(z); nf = np.max(np.abs(f))
        if np.max(np.abs(resid(z))) < 1e-14: return z, True
        J = np.zeros((2*m, 2*m), dtype=complex)
        for col in range(2*m):
            h = 1e-7*(1 + abs(z[col])); z2 = z.copy(); z2[col] += h
            J[:, col] = (dresid(z2) - f)/h
        try: dz = np.linalg.lstsq(J, -f, rcond=None)[0]
        except np.linalg.LinAlgError: return z, False
        lam = 1.0
        while lam > 1e-6:
            z2 = z + lam*dz
            if REAL: z2 = z2.real.astype(complex)
            if np.max(np.abs(dresid(z2))) < nf: break
            lam /= 2
        z = z2
        if np.max(np.abs(z)) > 1e5: return z, False
    return z, np.max(np.abs(resid(z))) < 1e-12
nsol = 0
for st in range(N):
    if REAL:
        V0 = rng.uniform(1/32, 3, m); U0 = 0.5 + (np.sqrt(8*V0) - 0.5 - 0.5)*rng.uniform(-1, 1, m)
        z0 = np.concatenate([U0, V0]).astype(complex)
    else:
        z0 = np.concatenate([rng.normal(0.5, 1.5, m) + 1j*rng.normal(0, 1, m), rng.uniform(0.05, 3, m)*np.exp(1j*rng.uniform(-1.5, 1.5, m))])
    z, ok = newton(z0)
    if not ok: continue
    U, W = z[:m], z[m:]
    if np.min(np.abs(W)) < 1e-3: continue
    D = np.abs(U[:, None] - U[None, :]) + np.abs(W[:, None] - W[None, :]) + np.eye(m)*1e9
    if D.min() < 1e-4: continue
    G = Gser(U, W); V = G[:k][::-1]    # V(t) = sum G_n t^(k-1-n); print low->high in t
    Vm1 = np.polyval(V[::-1], -1.0); dVm1 = np.polyval(np.polyder(V[::-1]), -1.0)
    cert = abs(Vm1 - 2.0**(1-k)) + abs(dVm1 - (1-k)*2.0**(1-k))
    if cert > 1e-8:
        print("  spurious st=%d cert=%.2e maxz=%.2e z=%s" % (st, cert, np.max(np.abs(z)), repr(list(z))), flush=True); continue
    nsol += 1
    isreal = np.max(np.abs(z.imag)) < 1e-8
    print("SOL st=%d real=%d V(t^0..t^%d)=%s S=%s T=%s plusres=%.2e" % (st, isreal, k-1, np.array2string(V, precision=5, max_line_width=500, suppress_small=True), np.round(np.sum(W), 6), np.round(np.sum(W*U), 6), np.max(np.abs(presid(z)))), flush=True)
print("m=%d lin=%d: %d nondegenerate solutions from %d starts" % (m, lin, nsol, N))
