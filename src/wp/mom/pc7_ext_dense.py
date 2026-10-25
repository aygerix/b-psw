# dense scan: is there an admissible c with ReX(th-) > 0 and ReY(th+) < 0 (bottom member condition)?
import numpy as np
def funcs(r1, r2, c):
    e1, e2 = r1 + r2, r1 * r2
    Dl = (1 - r1**2) * (1 - r2**2)
    lam, mu = 2 * e1 / Dl, 2 * (1 - e1**2 + e2) / Dl
    sm, qm = (c * e1 - lam) / c, (c * e2 + mu) / c
    sp, qp = -lam / c, (mu + lam * e1 - c) / c
    def two_re_X(s, q): return lam * (s**3 - 3 * q * s) + mu * (s**2 - 2 * q) - lam * s - 2 * mu
    RX = two_re_X(sm, qm); RY = two_re_X(sp, qp) + 4
    dm = sm**2 - 4 * qm; dp = sp**2 - 4 * qp
    return RX, RY, dm, dp
n = 0; bad = []
g = np.linspace(-0.999, 0.999, 401)
cs = np.exp(np.linspace(-12, 12, 6001))
worst = None
for i, r1 in enumerate(g):
    for r2 in g[i + 1:]:
        RX, RY, dm, dp = funcs(r1, r2, cs)
        adm = (dm < 0) & (dp < 0)
        if not adm.any(): continue
        n += 1
        m = adm & (RX > 0) & (RY < 0)
        if m.any(): bad.append((r1, r2, cs[m].min(), cs[m].max()))
        # margin: min over admissible of max(-RX, RY) scaled
        sc = np.maximum(-RX / (abs(RX) + abs(RY) + 1e-300), RY / (abs(RX) + abs(RY) + 1e-300))[adm]
        v = sc.min()
        if worst is None or v < worst[0]: worst = (v, r1, r2, cs[adm][sc.argmin()])
print("shapes with admissible c:", n, " bad:", len(bad))
print(bad[:20])
print("worst normalized margin", worst)
