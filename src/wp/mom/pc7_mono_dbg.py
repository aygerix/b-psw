import numpy as np
from pc7_mono_ext_lib import rho, shape
r1, r2 = -0.3, 0.6
A, B = 1 + r1 * r2, r1 + r2; Dl = A * A - B * B
e1, e2, lam, mu = shape(r1, r2)
print("Dl(shape) =", (1 - r1**2) * (1 - r2**2), " A^2-B^2 =", Dl)
lo, hi = A - np.sqrt(Dl), A + np.sqrt(Dl) - B * B / 2
print("I_x =", lo, hi)
for x in np.linspace(lo, hi, 5)[1:-1]:
    w = x * Dl / B**2; c = 1 / w
    dm = (lam - c * e1)**2 - 4 * c * (c * e2 + mu); dp = lam**2 - 4 * c * (mu + lam * e1 - c)
    psi = 2 * x**3 - (4 * A - B * B) * x**2 + (A * A - A * B * B + B * B) * x; tau = B * B * Dl / 2
    print(x, c, "discs", dm, dp, rho(c, r1, r2), "phi", psi / tau, " gp-form rho-:", (Dl**4 * c**3 + 2*Dl**2*(A*A - A*B*B + B*B)*c**2 + 2*B*B*Dl*(B*B - 4*A)*c + 4*B**4) / (Dl**4 * c**4))
