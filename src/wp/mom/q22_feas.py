"""Numerical feasibility of the 8 boundary equations (Theorem 12 setting) for (2,2) quartic pencils, k members, and of the
FULL Williams system. Unknowns: r1, r2 in (-1, 1) (via tanh), c_j > 0 (via exp). Boundary residuals: remainders of
V - 2^(1-k) G^(k-1), V' - ..., mod H, and W - 2^(1-k) H^(k-1), W' - ... mod G, with V = quotient((F - 1), G prod K-),
W = quotient((F + 1), H prod K+). Full residuals: remainders of (F - 1)/G mod K-_j and (F + 1)/H mod K+_j.
Also reports reality: K-+_j without real roots (needed), members distinct.
Usage: python3 q22_feas.py k mode(bnd|full) starts seed"""
import sys, warnings, numpy as np
warnings.filterwarnings('ignore')
from numpy.polynomial import polynomial as P
from scipy.optimize import least_squares
k = int(sys.argv[1]); mode = sys.argv[2]; starts = int(sys.argv[3]); rng = np.random.default_rng(int(sys.argv[4]))
def pm(a, b): return P.polymul(a, b)
def pdivmod(a, b):
    q, r = P.polydiv(a, b); return q, r
def setup(z):
    z = np.clip(z, -8, 8)
    r1, r2 = np.tanh(z[0]), np.tanh(z[1]); cs = np.exp(z[2:2 + k])
    G = np.array([-1.0, 0, 1]); H = pm([-r1, 1], [-r2, 1])
    Dl = (1 - r1**2) * (1 - r2**2); e1 = r1 + r2; e2 = r1 * r2
    lam = 2 * e1 / Dl; mu = 2 * (1 - e1**2 + e2) / Dl
    f0 = P.polyadd([1.0], pm(G, [mu, lam]))
    GH = pm(G, H)
    fs = [P.polyadd(f0, c * GH) for c in cs]
    return r1, r2, cs, G, H, fs
def resid(z):
    r1, r2, cs, G, H, fs = setup(z)
    F = np.array([1.0])
    for f in fs: F = pm(F, f)
    Km = [pdivmod(P.polysub(f, [1.0]), G)[0] for f in fs]
    Kp = [pdivmod(P.polyadd(f, [1.0]), H)[0] for f in fs]
    out = []
    if mode == 'bnd':
        LM = G.copy(); LP = H.copy()
        for K in Km: LM = pm(LM, K)
        for K in Kp: LP = pm(LP, K)
        V = pdivmod(P.polysub(F, [1.0]), LM)[0]; W = pdivmod(P.polyadd(F, [1.0]), LP)[0]
        c0 = 2.0 ** (1 - k); Gk = P.polypow(G, k - 1); Hk = P.polypow(H, k - 1)
        for A, B, M in [(V, c0 * Gk, H), (P.polyder(V), c0 * P.polyder(Gk), H), (W, c0 * Hk, G), (P.polyder(W), c0 * P.polyder(Hk), G)]:
            rr = pdivmod(P.polysub(A, B), M)[1]; rr = np.pad(rr, (0, 2 - len(rr)))
            out.extend(rr[:2])
    else:
        Fm = pdivmod(P.polysub(F, [1.0]), G)[0]; Fp = pdivmod(P.polyadd(F, [1.0]), H)[0]
        for K in Km:
            rr = pdivmod(Fm, K)[1]; rr = np.pad(rr, (0, 2 - len(rr))); out.extend(rr[:2] / max(1.0, np.abs(Fm).max()) )
        for K in Kp:
            rr = pdivmod(Fp, K)[1]; rr = np.pad(rr, (0, 2 - len(rr))); out.extend(rr[:2] / max(1.0, np.abs(Fp).max()))
    # distinctness penalty (soft): members must not coincide
    lc = np.log(cs); d = np.abs(lc[:, None] - lc[None, :])[np.triu_indices(k, 1)]
    out.extend(1e-3 / (d + 1e-3) - 0.0)  # pushes apart only when very close
    return np.array(out)
def reality(z):
    r1, r2, cs, G, H, fs = setup(z)
    bad = 0
    for f in fs:
        Km = pdivmod(P.polysub(f, [1.0]), G)[0]; Kp = pdivmod(P.polyadd(f, [1.0]), H)[0]
        for K in (Km, Kp):
            c, b, a = (list(K) + [0, 0, 0])[:3]
            if b * b - 4 * a * c >= 0: bad += 1
    return bad
best = []
for s in range(starts):
    z0 = np.concatenate([rng.uniform(-2, 2, 2), rng.uniform(-3, 3, k)])
    try:
        sol = least_squares(resid, z0, method='lm', max_nfev=600, xtol=1e-15, ftol=1e-15)
    except Exception as e:
        continue
    rv = resid(sol.x); core = rv[:-(k * (k - 1) // 2)]
    best.append((np.linalg.norm(core), sol.x))
    if s % 10 == 9: print(f'  start {s+1}: best so far {min(b[0] for b in best):.3e}', flush=True)
best.sort(key=lambda x: x[0])
print(f"k = {k}, mode = {mode}, starts = {starts}: residual quantiles", np.quantile([b[0] for b in best], [0, 0.1, 0.5, 0.9]) if best else None)
for nr, z in best[:5]:
    r1, r2, cs, *_ = setup(z)
    print(f"  res {nr:.3e}  r1 {r1:.6f} r2 {r2:.6f}  c {np.round(np.sort(cs), 5)}  reality-violations {reality(z)}  min|logc gap| {np.min(np.abs(np.subtract.outer(np.log(cs), np.log(cs)))[np.triu_indices(k,1)]):.2e}")
