#!/usr/bin/env python3
"""SUPERSEDED by ns_run3.py.  This v1 driver uses the per-member remainder residual only; that residual has
spurious zeros (diagonal, coincident special points, a -> infinity) and a flat landscape when |P_i| >> 1.
ns line: multistart batched Levenberg-Marquardt search for (real / complex) cubic Williams families.

Usage: ns_run.py k lin field sides nstarts seed [sigma] [batch] [iters]
  field: 'real'  -> exact parametrisation of the reality domain (g_i, h_i without real roots), 2m unknowns
         'free'  -> real (a, c) unrestricted, 2m unknowns
         'cx'    -> complex (a, c), 4m real unknowns
  sides: '-+' (full family, 4m equations) or '-' (the "-" half, 2m equations) or '+'
  sigma: deflation strength for coincident members (residual multiplied by 1 + sigma sum 1/d_ij^2), default 1e-3
Stage 1: deflated LM from random starts.  Stage 2: undeflated LM polish of every run (from the stage-1 endpoint),
to see whether the raw residual goes to 0 while members stay apart.  Output: summary to stdout, all endpoints to
ns_out/<tag>.npz.
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ns_core as nc


def make(k, lin, field, sides, sigma):
    m = k - lin

    def ac_of(x):
        if field == 'real':
            u, v = nc.real_unpack(x)
            a, c = nc.uv_to_ac(u, v)
            return a + 0j, c + 0j
        if field == 'free':
            return nc.free_unpack(x)
        if field == 'box':   # real, 0 < a < 40, c free
            m2 = x.shape[1] // 2
            return (40.0 / (1 + np.exp(-np.clip(x[:, :m2], -40, 40)))) + 0j, x[:, m2:] + 0j
        return nc.cx_unpack(x)

    def raw(x):
        a, c = ac_of(x)
        with np.errstate(all='ignore'):
            R = nc.full_resid(a, c, lin, sides)
        if field == 'cx':
            Rr = np.concatenate([R.real, R.imag], axis=1)
        else:
            Rr = R.real
        Rr = np.where(np.isfinite(Rr), Rr, 1e3)
        return Rr, a, c

    def fun(x):
        Rr, a, c = raw(x)
        if sigma > 0:
            D = nc.deflation(a, c, lin, sigma, sides=sides)
            Rr = Rr * D[:, None]
        return Rr

    return m, ac_of, raw, fun


def starts(rng, B, m, field):
    if field == 'real':
        a = np.exp(rng.uniform(np.log(0.05), np.log(31.5), (B, m)))
        v = 1 / a
        r = np.sqrt(8 * v) - 0.5
        return np.concatenate([np.log(r), rng.normal(0, 1.0, (B, m))], axis=1)
    if field == 'box':
        a = np.exp(rng.uniform(np.log(0.05), np.log(30), (B, m)))
        c = rng.normal(0, 1, (B, m)) * np.sqrt(8 * a) + 0.5 * a
        return np.concatenate([np.log(a / (40 - a)), c], axis=1)
    if field == 'free':
        a = np.exp(rng.uniform(np.log(0.05), np.log(30), (B, m))) * rng.choice([1, 1, 1, -1], (B, m))
        c = rng.normal(0, 1, (B, m)) * np.sqrt(8 * np.abs(a))
        return np.concatenate([a, c], axis=1)
    # complex: real-domain-like start plus complex noise of random size
    a = np.exp(rng.uniform(np.log(0.05), np.log(30), (B, m)))
    c = rng.uniform(-1, 1, (B, m)) * np.sqrt(8 * a) * 0.7 + 0.5 * a
    sa = np.exp(rng.uniform(np.log(0.01), np.log(3), (B, 1)))
    ai = rng.normal(0, 1, (B, m)) * sa * a; ci = rng.normal(0, 1, (B, m)) * sa * (np.abs(c) + 1)
    return np.concatenate([a, ai, c, ci], axis=1)


def main():
    k, lin, field, sides = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
    ns, seed = int(sys.argv[5]), int(sys.argv[6])
    sigma = float(sys.argv[7]) if len(sys.argv) > 7 else 1e-3
    B = int(sys.argv[8]) if len(sys.argv) > 8 else 500
    iters = int(sys.argv[9]) if len(sys.argv) > 9 else 300
    tag = f"k{k}_l{lin}_{field}_{'m' if sides == '-' else ('p' if sides == '+' else 'mp')}_s{seed}_sg{sigma:g}"
    m, ac_of, raw, fun = make(k, lin, field, sides, sigma)
    _, _, _, fun0 = make(k, lin, field, sides, 0.0)
    rng = np.random.default_rng(seed)
    X1, X2, C1 = [], [], []
    t0 = time.time()
    done = 0
    while done < ns:
        b = min(B, ns - done)
        x0 = starts(rng, b, m, field)
        x1, c1 = nc.batched_lm(fun, x0, iters=iters)
        x2, c2 = nc.batched_lm(fun0, x1, iters=iters)
        X1.append(x1); X2.append(x2); C1.append(c1)
        done += b
        print(f"[{tag}] {done}/{ns} starts, {time.time() - t0:.0f}s", flush=True)
    X1 = np.concatenate(X1); X2 = np.concatenate(X2)
    out = {}
    for name, X in (('defl', X1), ('polish', X2)):
        Rr, a, c = raw(X)
        res = np.abs(Rr).max(axis=1)
        sep = nc.pair_sep(a, c, lin)
        mina = np.abs(a).min(axis=1)
        mg = nc.margins(a, c) if field != 'cx' else np.full(len(X), np.nan)
        with np.errstate(all='ignore'):
            gap = nc.root_gap_terms(a, c, sides)[1]
        imag = (np.abs(a.imag).max(axis=1) + np.abs(c.imag).max(axis=1)) if field == 'cx' else np.zeros(len(X))
        out[name] = (res, sep, mina, mg, imag, a, c)
        np.savez_compressed(f"ns_out/{tag}_{name}.npz", x=X, res=res, sep=sep, mina=mina, mg=mg, imag=imag,
                            gap=gap, a=a, c=c)
        q = np.quantile(res, [0, 0.01, 0.1, 0.5, 0.9])
        print(f"== {tag} [{name}] max-abs normalised residual quantiles (0,1,10,50,90%): "
              + " ".join(f"{z:.3e}" for z in q))
        for dl in (1e-6, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.3, 1.0):
            ok = (sep >= dl) & (gap >= 1e-3)
            if ok.any():
                i = np.argmin(np.where(ok, res, np.inf))
                print(f"   sep >= {dl:<6g} & rootgap >= 1e-3: n = {ok.sum():6d}, min residual {res[i]:.3e} (sep {sep[i]:.3e}, "
                      f"gap {gap[i]:.3e}, min|a| {mina[i]:.3e}, margin {mg[i]:.3e}, imag {imag[i]:.2e})")
        for thr in (1e-12, 1e-9, 1e-6, 1e-3):
            z = res < thr
            print(f"   residual < {thr:g}: {z.sum()} runs; of those sep >= 1e-3: {(z & (sep >= 1e-3)).sum()}"
                  f"; and rootgap >= 1e-3: {(z & (sep >= 1e-3) & (gap >= 1e-3)).sum()}")
        good = (sep >= 1e-2) & (gap >= 1e-2)
        order = np.argsort(np.where(good, res, np.inf))[:3]
        for i in order:
            if not np.isfinite(res[i]) or not good[i]:
                continue
            print(f"   best: res {res[i]:.3e} sep {sep[i]:.3e} gap {gap[i]:.3e} margin {mg[i]:.3e}\n      a = {np.array2string(a[i], precision=8, max_line_width=300)}"
                  f"\n      c = {np.array2string(c[i], precision=8, max_line_width=300)}")
    print(f"total time {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
