"""Exact-form extreme-member test (no truncation in c). A = 1 + r1 r2, B = r1 + r2, Dl = A^2 - B^2,
x = B^2 w / Dl (w = 1/c), reality interval I = (A - sqrt(Dl), A + sqrt(Dl) - B^2/2),
psi(x) = 2x^3 - (4A - B^2)x^2 + (A^2 - A B^2 + B^2)x, tau = B^2 Dl / 2.
Max member needs psi(x_b) > tau, min member needs psi(x_a) < -tau, with x_b < x_a in I.  Feasible iff inf S+ < sup S-.
Also cross-checks rho+- from pc7_mono_ext_lib at a few c."""
import numpy as np, sys
from pc7_mono_ext_lib import rho
def test(r1, r2, n=20001):
    A, B = 1 + r1 * r2, r1 + r2; Dl = A * A - B * B; s = np.sqrt(Dl)
    lo, hi = A - s, A + s - B * B / 2
    if lo >= hi: return 'empty', None
    x = np.linspace(lo, hi, n)[1:-1]
    psi = 2 * x**3 - (4 * A - B * B) * x**2 + (A * A - A * B * B + B * B) * x; tau = B * B * Dl / 2
    Sp, Sm = x[psi > tau], x[psi < -tau]
    feas = len(Sp) > 0 and len(Sm) > 0 and Sp.min() < Sm.max()
    return ('feasible' if feas else 'excluded'), (lo, hi)
# cross-check against direct rho at random shape / c in the reality interval
rng = np.random.default_rng(3); mism = 0
for _ in range(200):
    r1, r2 = np.sort(rng.uniform(-1, 1, 2)); A, B = 1 + r1 * r2, r1 + r2; Dl = A * A - B * B
    if B * B < 1e-6: continue
    lo, hi = A - np.sqrt(Dl), A + np.sqrt(Dl) - B * B / 2
    if lo >= hi: continue
    x = rng.uniform(lo, hi); w = x * Dl / B**2; c = 1 / w
    ok, rm, rp = rho(c * 1.0, r1, r2)   # NOTE: c here is the pencil coefficient; Dl in w relation is Dl of shape
    psi = 2 * x**3 - (4 * A - B * B) * x**2 + (A * A - A * B * B + B * B) * x; tau = B * B * Dl / 2
    phi = psi / tau
    if not ok or np.sign(rm) != np.sign(1 + phi) or np.sign(rp) != np.sign(1 - phi): mism += 1
print("cross-check mismatches:", mism)
cnt = {}; feas = []
for r1 in np.linspace(-0.99, 0.99, 100):
    for r2 in np.linspace(-0.99, 0.99, 100):
        if r2 <= r1 + 1e-9: continue
        st, I = test(r1, r2, 4001); cnt[st] = cnt.get(st, 0) + 1
        if st == 'feasible': feas.append((r1, r2))
print(cnt)
F = np.array(feas)
if len(F): print("feasible |e1| range:", np.abs(F.sum(1)).min(), np.abs(F.sum(1)).max(), " sample:", F[:: max(1, len(F)//10)].round(3).tolist())
