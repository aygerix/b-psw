"""Extreme-member sign test for real (2,2) pencils, k members, c_j > 0 (forced by reality).
For the member with the SMALLEST c: V-half at the roots z of K-(c) = c H + L forces Re(GH(z)) < 0, W-half at the roots of
K+(c) = c G + L' forces Re(GH(z')) > 0. For the LARGEST c the signs flip. rho-(c) = Re GH at roots of K-, rho+(c) = Re GH at roots of K+.
Scan: reality interval I_c (both K-+ without real roots); sign patterns of rho-+ on I_c."""
import numpy as np, sys
def shape(r1, r2):
    e1, e2 = r1 + r2, r1 * r2; Dl = (1 - r1**2) * (1 - r2**2)
    return e1, e2, 2 * e1 / Dl, 2 * (1 - e1**2 + e2) / Dl
def GH(z, e1, e2): return (z * z - 1) * (z * z - e1 * z + e2)
def rho(c, r1, r2):
    e1, e2, lam, mu = shape(r1, r2)
    Km = np.roots([c, lam - c * e1, c * e2 + mu]); Kp = np.roots([c, lam, mu + lam * e1 - c])
    ok = np.all(np.abs(Km.imag) > 1e-12) and np.all(np.abs(Kp.imag) > 1e-12)
    return ok, GH(Km[0], e1, e2).real, GH(Kp[0], e1, e2).real
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
stats = {}; bad = []
for it in range(int(sys.argv[2]) if len(sys.argv) > 2 else 400):
    r1, r2 = np.sort(rng.uniform(-1, 1, 2))
    cs = np.exp(np.linspace(-6, 6, 1200)); data = [(c,) + rho(c, r1, r2) for c in cs]
    good = [(c, a, b) for c, ok, a, b in data if ok]
    if not good: stats['empty'] = stats.get('empty', 0) + 1; continue
    C = np.array(good)
    # is there c1 < c2 in I_c with rho-(c1) < 0 < rho-(c2) and rho+(c1) > 0 > rho+(c2)?
    neg_m = C[:, 1] < 0; pos_m = C[:, 1] > 0; pos_p = C[:, 2] > 0; neg_p = C[:, 2] < 0
    A = neg_m & pos_p    # admissible for the minimal member
    B = pos_m & neg_p    # admissible for the maximal member
    feas = A.any() and B.any() and C[A, 0].min() < C[B, 0].max()
    key = 'feasible' if feas else 'excluded'
    stats[key] = stats.get(key, 0) + 1
    if feas and len(bad) < 10: bad.append((r1, r2, C[A, 0].min(), C[B, 0].max()))
print(stats)
for b in bad: print("  feasible shape r = (%.4f, %.4f): min-admissible c from %.4g, max-admissible c up to %.4g" % b)
