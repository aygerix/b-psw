# pc7_adv_num.py -- discriminating numerical test of the Pencil Lemma (adversarial verifier).
# General complex pencil, degree d, k members: P monic (complex), f0 deg < d (complex), real c_1 < ... < c_k.
# Impose the Williams conditions ONLY at member 1 (the smallest c):
#   at roots th of f_1 - 1:  prod_{i>1} (1 + a_i P(th)) = 1
#   at roots th' of f_1 + 1: (-1)^k prod_{i>1} (1 - a_i P(th')) = -1
# a_i = c_i - c_1 > 0. The lemma says: no solution with c_1 > 0. Positive control: c_1 < 0 should be solvable
# (system is underdetermined), so the test can tell the two cases apart.
import numpy as np, sys
from scipy.optimize import least_squares
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 1)

def resid(z, d, k, sgn):
    P = np.concatenate([[1.0 + 0j], z[0:d] + 1j * z[d:2 * d]])           # monic, highest first
    f0 = np.concatenate([[0j], z[2 * d:3 * d] + 1j * z[3 * d:4 * d]])     # deg < d
    c1 = sgn * np.exp(z[4 * d]); a = np.exp(z[4 * d + 1:4 * d + k])
    f1 = f0 + c1 * P
    fm = f1.copy(); fm[-1] -= 1; fp = f1.copy(); fp[-1] += 1
    out = []
    for th in np.roots(fm):
        w = np.polyval(P, th); out.append(np.prod(1 + a * w) - 1)
    for th in np.roots(fp):
        w = np.polyval(P, th); out.append((-1) ** k * np.prod(1 - a * w) + 1)
    out = np.array(out)
    return np.concatenate([out.real, out.imag])

for d, k in [(1, 3), (2, 3), (2, 5), (3, 3), (4, 3), (4, 7)]:
    for sgn in (+1, -1):
        best = np.inf; nsol = 0; bestz = None
        for s in range(100):
            z0 = rng.normal(size=4 * d + k) * 0.8
            try:
                r = least_squares(resid, z0, args=(d, k, sgn), method='trf', max_nfev=3000, xtol=1e-15, ftol=1e-15)
            except Exception:
                continue
            v = np.sqrt(np.sum(r.fun ** 2))
            # reject runaway: exp params too large (c or a -> 0 / inf degenerate limits)
            if np.max(np.abs(r.x[4 * d:])) > 12: continue
            if v < best: best, bestz = v, r.x
            if v < 1e-10: nsol += 1
        info = ""
        if bestz is not None:
            info = "c1=%.4g a=%s" % (sgn * np.exp(bestz[4 * d]), np.array2string(np.exp(bestz[4 * d + 1:]), precision=3))
        print("d=%d k=%d sign(c1)=%+d  best residual=%.3e  #solutions(<1e-10)=%d  %s" % (d, k, sgn, best, nsol, info), flush=True)
