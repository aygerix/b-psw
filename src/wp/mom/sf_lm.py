# sf_lm.py -- discriminating experiment for the sign-free extension (Lemma 3, Theorem 10, Theorem 11):
# solve the FULL Williams system for REAL members with a PRESCRIBED leading-sign pattern, and as positive control over
# C (where zeta-families {zeta^j q} exist for odd k).
# Williams <=> for every member i and e = +-1: P_i(rho) = 1 at every root rho of f_i - e (P_i = prod_{l != i} f_l),
# with multiplicity; the solver uses the root values (generic simple roots) and the final point is re-checked with the
# exact remainder test rem(F - e, f_i - e) (polydiv).
# Gauge: f_1 = s_1 t^d + 0 t^(d-1) + ... (scale and translation of t).  Leading coefficients of the other members
# s_i exp(alpha_i) (real) or exp(alpha_i + i beta_i) (complex), alpha clipped to [-12, 12].
# Deflation of coincident members (optional 'defl'): residual * (1 + 1/min_{i<j} |f_i - f_j|^2) (same degree pairs).
# Usage: python3 sf_lm.py MODE DEGS NNEG STARTS SEED [defl] [NFEV_PER_UNKNOWN]
import sys, numpy as np
from scipy.optimize import least_squares
np.seterr(all="ignore")
mode, degs, nneg, starts, seed = sys.argv[1], [int(z) for z in sys.argv[2].split(",")], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
DEFL = len(sys.argv) > 6 and sys.argv[6] == "defl"
NFEV = int(sys.argv[7]) if len(sys.argv) > 7 else 200
k = len(degs); rng = np.random.default_rng(seed)
signs = np.array([1]*(k - nneg) + [-1]*nneg)
cplx = (mode == "cplx")
def nfree():
    n = max(degs[0] - 1, 0)
    for d in degs[1:]: n += 1 + d
    return n * (2 if cplx else 1)
def unpack(x):
    z = (x[:len(x)//2] + 1j*x[len(x)//2:]) if cplx else x
    mem = []; p = 0; dt = complex if cplx else float
    d = degs[0]; c = np.zeros(d + 1, dtype=dt); c[0] = 1.0 if cplx else signs[0]
    if d >= 2: c[2:] = z[p:p + d - 1]; p += d - 1
    mem.append(c)
    for i in range(1, k):
        d = degs[i]; c = np.zeros(d + 1, dtype=dt)
        a = np.clip(z[p].real, -12, 12) + (1j*z[p].imag if cplx else 0)
        c[0] = np.exp(a) if cplx else signs[i]*np.exp(a); p += 1
        c[1:] = z[p:p + d]; p += d
        mem.append(c)
    return mem
def roots_of(c):
    d = len(c) - 1
    if d == 1: return np.array([-c[1]/c[0]], dtype=complex)
    if d == 2:
        A, B, C = complex(c[0]), complex(c[1]), complex(c[2]); s = np.sqrt(B*B - 4*A*C)
        q = -0.5*(B + (s if (B.conjugate()*s).real >= 0 else -s))
        return np.array([q/A, C/q if q != 0 else -B/A])
    return np.roots(c).astype(complex)
def resid_root(mem):
    out = []
    for i, c in enumerate(mem):
        for e in (1.0, -1.0):
            D = c.astype(complex).copy(); D[-1] -= e
            rho = roots_of(D)
            P = np.ones(len(rho), dtype=complex)
            for l, cl in enumerate(mem):
                if l != i: P = P * np.polyval(cl, rho)
            out.append(P - 1.0)
    out = np.concatenate(out)
    return np.concatenate([out.real, out.imag])
def resid_rem(mem):   # exact remainder test, for the final check
    F = np.array([1.0 + 0j])
    for c in mem: F = np.polymul(F, c)
    out = []
    for c in mem:
        d = len(c) - 1
        for e in (1.0, -1.0):
            D = c.astype(complex).copy(); D[-1] -= e; D = D / D[0]
            Fe = F.copy(); Fe[-1] -= e
            _, r = np.polydiv(Fe, D)
            r = np.concatenate([np.zeros(d - len(r)), r]) if len(r) < d else r[-d:]
            out.append(r)
    out = np.concatenate(out); return np.concatenate([out.real, out.imag])
def dmin2(mem):
    d = np.inf
    for i in range(k):
        for j in range(i + 1, k):
            if len(mem[i]) == len(mem[j]): d = min(d, float(np.sum(np.abs(mem[i] - mem[j])**2)))
    return d
def fun(x):
    mem = unpack(x); r = resid_root(mem)
    if DEFL:
        d = dmin2(mem)
        if np.isfinite(d): r = r * (1.0 + 1.0/max(d, 1e-300))
    r[~np.isfinite(r)] = 1e8
    return r
n = nfree()
stats = dict(zero_distinct=0, zero_coalesced=0, zero_degenerate=0, nonzero=0)
best = (np.inf, None)
for s in range(starts):
    x0 = rng.normal(size=n) * rng.choice([0.3, 1.0, 2.0])
    try:
        r = least_squares(fun, x0, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=NFEV * (n + 1))
    except Exception as ex:
        stats["nonzero"] += 1; print("   start", s, "error", ex, flush=True); continue
    mem = unpack(r.x); rr = resid_root(mem); c2 = float(rr @ rr); rm = resid_rem(mem); c2rem = float(rm @ rm)
    leads = np.array([abs(m[0]) for m in mem]); big = max(np.max(np.abs(m)) for m in mem); dmin = np.sqrt(dmin2(mem))
    print(f"   start {s}: res2(root) = {c2:.3e}  res2(rem) = {c2rem:.3e}  dmin = {dmin:.3e}  |lead| in [{leads.min():.2e}, {leads.max():.2e}]  max|coef| = {big:.2e}  nfev = {r.nfev}", flush=True)
    if c2 < best[0]: best = (c2, (dmin, leads.min(), leads.max(), big, [np.round(m, 5) for m in mem], c2rem))
    if c2 < 1e-20 and c2rem < 1e-16:
        if leads.min() < 1e-5 or big > 1e5: stats["zero_degenerate"] += 1
        elif dmin < 1e-6: stats["zero_coalesced"] += 1
        else: stats["zero_distinct"] += 1
    else:
        stats["nonzero"] += 1
print(f"SUMMARY mode={mode} defl={DEFL} degs={degs} nneg={nneg} unknowns={n} starts={starts}: {stats}")
if best[1] is not None:
    print(f"   BEST res2(root) = {best[0]:.3e} res2(rem) = {best[1][5]:.3e}; min pair dist = {best[1][0]:.3e}, |lead| in [{best[1][1]:.3e}, {best[1][2]:.3e}], max |coef| = {best[1][3]:.3e}")
    if best[0] < 1e-20: print("   BEST members:", best[1][4])
sys.stdout.flush()
