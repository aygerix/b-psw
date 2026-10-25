#!/usr/bin/env python3
"""m9v: numerical cross-check of m9_hier.taylor_resid r_2 (normalized by c0 = 2^(1-k)) against
[e_2(s) - (2^k V''(-1) - 2(k-1)(k-2))]/4, plus r_0, r_1 (check mode of the claim's evidence)."""
import sys
import numpy as np
sys.path.insert(0, '..')
from m9_hier import taylor_resid, Vcoef
rng = np.random.default_rng(7)
for k, lin in [(9, 0), (9, 1), (11, 0), (11, 1)]:
    m = k - lin
    worst = 0.0
    for _ in range(20):
        u = 0.5 + rng.normal(0, .3, m)
        v = (np.abs(u - .5) + .5) ** 2 / 8 + rng.uniform(0.01, 1, m)
        r, _ = taylor_resid(u, v, k, lin, 3)
        c0 = 2.0 ** (1 - k)
        Vc = Vcoef(u, v, k, lin)
        P = np.polynomial.Polynomial(Vc)
        V0, V1, V2 = P(-1.0), P.deriv(1)(-1.0), P.deriv(2)(-1.0)
        s = 2 + (1 - u) / v
        if lin:
            s = np.append(s, 2.0)
        e2 = (s.sum() ** 2 - (s ** 2).sum()) / 2
        pred = np.array([c0 - V0, (1 - k) * c0 - V1, (e2 - (2 ** k * V2 - 2 * (k - 1) * (k - 2))) / 2 ** (k + 1)])
        rel = np.abs(r - pred) / (np.abs(pred) + 1e-300)
        worst = max(worst, rel.max())
        # normalized form of the claim: r_2 / c0 = [...] / 4
        assert abs(r[2] / c0 - (e2 - (2 ** k * V2 - 2 * (k - 1) * (k - 2))) / 4) <= 1e-9 * (abs(r[2] / c0) + 1)
    print(f"k={k} lin={lin}: worst relative discrepancy r0..r2 vs formulas {worst:.2e}")
