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
