# pc7_ext_num.py -- independent floating-point controls for the Pencil Lemma (numpy roots, no gp formulas reused)
import numpy as np
rng = np.random.default_rng(7)
P_ = np.polynomial.polynomial
# (A) root-sum identity sum_{f=1} P - sum_{f=-1} P = 2d/c for random f (deg d, lead c) and random monic P (deg d), complex too
worst = 0.0
for trial in range(20000):
    d = rng.integers(1, 9)
    c = rng.normal() * np.exp(rng.uniform(-2, 2))
    f = np.concatenate([[c], rng.normal(size=d)])            # highest degree first
    if rng.random() < 0.5: f = f + 1j * np.concatenate([[0], rng.normal(size=d)])
    P = np.concatenate([[1.0], rng.normal(size=d)])
    fm = f.copy(); fm[-1] -= 1; fp = f.copy(); fp[-1] += 1
    s = np.polyval(P, np.roots(fm)).sum() - np.polyval(P, np.roots(fp)).sum()
    err = abs(s - 2 * d / c) / max(1, abs(2 * d / c), np.abs(np.roots(fm)).max()**d)
    worst = max(worst, err)
print("(A) max normalized error of root-sum identity over 20000 random cases:", worst)
# (B) (2,2) pencil: random shape and c, ReP(th-) - ReP(th+) = 4/c from the actual quartic f_c
def pencil(r1, r2, c):
    e1, e2 = r1 + r2, r1 * r2; Dl = (1 - r1**2) * (1 - r2**2)
    lam, mu = 2 * e1 / Dl, 2 * (1 - e1**2 + e2) / Dl
    G = np.array([1.0, 0, -1]); H = np.array([1.0, -e1, e2]); L = np.array([lam, mu])
    f0 = np.polyadd([1.0], np.polymul(G, L)); P = np.polymul(G, H)
    return np.polyadd(f0, c * P), P, G, H
worst = 0.0; nadm = 0
for trial in range(20000):
    r1, r2 = np.sort(rng.uniform(-1, 1, 2)); c = np.exp(rng.uniform(-6, 6))
    f, P, G, H = pencil(r1, r2, c)
    fm = f.copy(); fm[-1] -= 1; fp = f.copy(); fp[-1] += 1
    rm, rp = np.roots(fm), np.roots(fp)
    # check the shared points are roots
    assert min(abs(rm - 1)) < 1e-6 and min(abs(rm + 1)) < 1e-6 and min(abs(rp - r1)) < 1e-5 * max(1, 1/(r2 - r1 + 1e-9))
    vm, vp = np.polyval(P, rm), np.polyval(P, rp)
    scale = max(8 / c, np.abs(vm).sum() + np.abs(vp).sum(), np.abs(np.concatenate([rm, rp])).max()**4)
    e = abs((vm.sum() - vp.sum()) - 8 / c) / scale
    if e > worst: worst, wcase = e, (r1, r2, c, vm.sum() - vp.sum(), 8 / c)
print("(B) (2,2) pencils: max scaled error of sum_- P - sum_+ P = 8/c over 20000 random (shape, c):", worst, "case", wcase)
# (C) half-plane lemma: nonzero roots u of prod_i (1 + a_i u) = 1, a_i > 0, have Re u < 0 ; of prod(1 - a_i u) = +-1 have Re u > 0
mx1 = -np.inf; mn2 = np.inf
for trial in range(20000):
    m = rng.integers(1, 12); a = np.exp(rng.uniform(-4, 4, m))
    poly = np.array([1.0])
    for ai in a: poly = np.polymul(poly, [ai, 1.0])
    q = poly.copy(); q[-1] -= 1; rts = np.roots(q); rts = rts[abs(rts) > 1e-9 * max(1, abs(rts).max())]
    if len(rts): mx1 = max(mx1, (rts.real / abs(rts)).max())
    poly2 = np.array([1.0])
    for ai in a: poly2 = np.polymul(poly2, [-ai, 1.0])
    for sgn in (1, -1):
        q2 = poly2.copy(); q2[-1] -= sgn; rts2 = np.roots(q2); rts2 = rts2[abs(rts2) > 1e-9 * max(1, abs(rts2).max())]
        if len(rts2): mn2 = min(mn2, (rts2.real / abs(rts2)).min())
print("(C) max Re(u)/|u| over nonzero roots of prod(1 + a u) = 1:", mx1, ";  min Re(u)/|u| for prod(1 - a u) = +-1:", mn2)
# (D) sign forcing in the (2,2) setting: for -1 < r1 < r2 < 1 and c < 0, K- = cH + L has a real root or f - 1 > 0 at r1 (f(r1) = -1 contradiction)
bad = 0
for trial in range(20000):
    r1, r2 = np.sort(rng.uniform(-1, 1, 2)); c = -np.exp(rng.uniform(-6, 6))
    f, P, G, H = pencil(r1, r2, c)
    e1, e2 = r1 + r2, r1 * r2; Dl = (1 - r1**2) * (1 - r2**2); lam, mu = 2 * e1 / Dl, 2 * (1 - e1**2 + e2) / Dl
    Km = np.polyadd(c * H, [lam, mu]); disc = Km[1]**2 - 4 * Km[0] * Km[2]
    if disc < 0: bad += 1     # would be a counterexample to the sign-forcing claim
print("(D) c < 0 samples with K-_c free of real roots (claim: 0):", bad)
