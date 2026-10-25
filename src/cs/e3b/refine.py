#!/usr/bin/env python3
"""src/cs/e3b/refine.py -- E3b: post-processing of stratum files from strat_lp.py.
For each kept solution (c, l, r): (1) detect Lucas points with a loose tolerance (1e-5), (2) re-impose ALL detected
Lucas points in the Lucas-point model and re-converge by LM from the solution (residual < 1e-12 required), (3) recompute
the gcd degrees at tolerance 1e-9, (4) local dimension by slicing: the largest d such that LM with d random affine
slices through a perturbation of the point re-converges (component dimension >= d), (5) cluster.
Usage: refine.py k infile [maxpts]"""
import sys, numpy as np
from collections import Counter
sys.path.insert(0, __file__.rsplit('/', 1)[0] if '/' in __file__ else '.')
from strat_lp import LP, lm_batch, to_clr
from strat_batch import gcd_degrees, nondeg

def lucas_table(clr, k):
    c = clr[0]; lam = clr[1:k]; rr = clr[k:]
    g = lambda i, t: c if i == 0 else lam[i - 1] * (t - rr[i - 1])
    tab = []
    for j in range(k):
        Ncoef = [c, 0, 2 - c] if j == 0 else [lam[j-1], -lam[j-1]*rr[j-1], -lam[j-1], lam[j-1]*rr[j-1] + 2]
        rts = np.roots(Ncoef); row = []
        for t in rts:
            val = np.prod([g(j, t) - 2 * g(i, t) for i in range(k) if i != j]) - g(j, t)**(k - 1)
            row.append((t, abs(val) / abs(g(j, t))**(k - 1)))
        tab.append(row)
    return tab

def build_z(model, clr, pts):
    """Lucas-point-model vector from (c, l, r) and Lucas points per member (lists)."""
    k = model.k; z = np.zeros(model.n, dtype=complex); z[:k - 1] = clr[k:]
    for j in range(k):
        kind, ix = model.lay[j]
        if kind == 'coef':
            z[ix] = clr[0] if j == 0 else clr[j]
        else:
            P = list(pts[j])
            if j == 0:
                z[ix[0]] = P[0]
            else:
                for t, idx in enumerate(ix):
                    z[idx] = P[t]
    return z

def slice_dim(model_args, z, k, rng, dmax=5):
    best = 0
    for d in range(1, dmax + 1):
        okc = 0
        for trial in range(3):
            mdl = LP(*model_args, nslice=d, rng=rng)
            zp = z + 1e-3 * (rng.normal(size=z.size) + 1j * rng.normal(size=z.size))
            mdl.s0 = (mdl.S @ zp)                     # slices through the perturbed point
            zz, nf = lm_batch(mdl, (z + 0j)[None].copy(), iters=300)
            if nf[0] < 1e-11 and np.linalg.norm(zz[0] - z) < 0.1:
                okc += 1
        if okc == 0:
            break
        best = d
    return best

def main():
    k = int(sys.argv[1]); fn = sys.argv[2]; maxpts = int(sys.argv[3]) if len(sys.argv) > 3 else 400
    rng = np.random.default_rng(7)
    rows = [eval(l) for l in open(fn)][:maxpts]
    out = []
    for d0, gd0, clr, pts in rows:
        clr = np.array(clr, dtype=complex)
        tab = lucas_table(clr, k)
        lp = [[t for (t, r) in row if r < 1e-5] for row in tab]
        m = [len(p) for p in lp]
        if m[0] == 1 and False:
            pass
        # member 1 with 2 points: the LP model uses th and -th
        model = LP(k, m, rng=rng)
        z = build_z(model, clr, lp)
        zz, nf = lm_batch(model, z[None].copy(), iters=300)
        if nf[0] > 1e-12:
            out.append(('norefine', m, None, clr)); continue
        clr2 = to_clr(model, zz[0])
        if not nondeg(clr2, k):
            out.append(('degenerate', m, None, clr2)); continue
        gd, _ = gcd_degrees(clr2, k, tol=1e-9)
        dim = slice_dim((k, m), zz[0], k, rng)
        out.append(('ok', gd, dim, clr2))
    cnt = Counter((s, tuple(g), d) for s, g, d, c in out)
    for key, v in sorted(cnt.items(), key=lambda t: (t[0][0], -sum(t[0][1]))):
        print(f"  {key[0]:10s} gcd {list(key[1])} sum {sum(key[1])}  dim {key[2]}: {v}")
    best = max([sum(g) for s, g, d, c in out if s == 'ok'], default=None)
    print("max refined sum m:", best)
    for s, g, d, c in out:
        if s == 'ok' and sum(g) == best:
            print("  rep", g, d, np.round(c, 6).tolist()); break

if __name__ == "__main__":
    main()
