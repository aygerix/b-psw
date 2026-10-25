#!/usr/bin/env python3
"""SUPERSEDED by ns_run3.py (only _lm_fixed_last is still imported from here).  The log residual helps the
global phase but evaluation at roots misses derivative conditions at coincidences.
ns line, driver v2.  Usage:
  ns_run2.py k lin field sides nstarts seed [sigma=1e-3] [batch=500] [it1=400] [it2=200] [planted=0]
field: real (exact reality-domain parametrisation) | cx (complex a, c) | box (real, 0 < a < 40, c free)
sides: -+ (full family) | - (half)
Stage 1: LOG residual times deflation (pairs of members, all special-point coincidences incl. 0, -1).
Stage 2: exact REMAINDER residual, no deflation (polish; shows whether the raw residual really goes to 0).
planted = 1: every start gets its own random target x* (tau = remainders of x*), a recall control: x* is an exact
zero of the planted system with the same structure.
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ns_core as nc
import ns_core2 as n2


def ac_of(x, field):
    if field == 'real':
        u, v = nc.real_unpack(x)
        a, c = nc.uv_to_ac(u, v)
        return a + 0j, c + 0j
    if field == 'box':
        m2 = x.shape[1] // 2
        return (40.0 / (1 + np.exp(-np.clip(x[:, :m2], -40, 40)))) + 0j, x[:, m2:] + 0j
    return nc.cx_unpack(x)


def starts(rng, B, m, field):
    if field == 'real':
        a = np.exp(rng.uniform(np.log(0.05), np.log(31.5), (B, m)))
        r = np.sqrt(8 / a) - 0.5
        return np.concatenate([np.log(r), rng.normal(0, 1.0, (B, m))], axis=1)
    if field == 'box':
        a = np.exp(rng.uniform(np.log(0.05), np.log(30), (B, m)))
        c = rng.uniform(-1, 1, (B, m)) * np.sqrt(8 * a) * 0.9 + 0.5 * a
        return np.concatenate([np.log(a / (40 - a)), c], axis=1)
    a = np.exp(rng.uniform(np.log(0.05), np.log(30), (B, m)))
    c = rng.uniform(-1, 1, (B, m)) * np.sqrt(8 * a) * 0.9 + 0.5 * a
    sa = np.exp(rng.uniform(np.log(0.01), np.log(2), (B, 1)))
    ai = rng.normal(0, 1, (B, m)) * sa * a; ci = rng.normal(0, 1, (B, m)) * sa * (np.abs(c) + 1)
    return np.concatenate([a, ai, c, ci], axis=1)


def real_rows(R, field):
    if field == 'cx':
        R = np.concatenate([R.real, R.imag], axis=1)
    else:
        R = R.real
    return np.where(np.isfinite(R), R, 1e3)


def build(k, lin, field, sides, sigma, tau_all):
    """tau_all: None or dict side -> (tau0, tau1) of shape (Btot, m); functions take (x, rows)"""
    def sub(rows):
        if tau_all is None:
            return None
        return {s: (tau_all[s][0][rows], tau_all[s][1][rows]) for s in sides}

    def flog(x, rows):
        a, c = ac_of(x, field)
        with np.errstate(all='ignore'):
            R = n2.log_resid(a, c, lin, sides, field, sub(rows))
            R = real_rows(R, field)
            if sigma > 0:
                D = n2.deflation(a, c, lin, sigma, sides)
                D = np.where(np.isfinite(D), D, 1e300)
                R = R * D[:, None]
        return R

    def frem(x, rows):
        a, c = ac_of(x, field)
        with np.errstate(all='ignore'):
            R = n2.rem_resid(a, c, lin, sides, sub(rows))
        return real_rows(R, field)
    return flog, frem


def lm_rows(f, x, rows, iters):
    """batched LM where f also needs the row indices (for per-row planted targets)"""
    # wrap: ns_core.batched_lm calls fun(x_subset) on subsets of rows; we track rows by appending an index column
    xi = np.concatenate([x, rows[:, None].astype(float)], axis=1)

    def g(z):
        r = z[:, -1].astype(int)
        return f(z[:, :-1], r)
    # the index column must not move: give it a zero Jacobian by evaluating with h -> exact same rows
    return _lm_fixed_last(g, xi, iters)


def _lm_fixed_last(fun, x0, iters):
    n = x0.shape[1] - 1
    x = x0.copy(); B = x.shape[0]
    R = fun(x); cost = np.einsum('bi,bi->b', R, R); lam = np.full(B, 1e-2)
    active = np.isfinite(cost); cost[~active] = np.inf
    for it in range(iters):
        ia = np.nonzero(active)[0]
        if len(ia) == 0:
            break
        xa = x[ia]; Ra = R[ia]
        J = np.empty((len(ia), Ra.shape[1], n))
        for col in range(n):
            hh = 1e-7 * (1 + np.abs(xa[:, col]))
            x2 = xa.copy(); x2[:, col] += hh
            J[:, :, col] = (fun(x2) - Ra) / hh[:, None]
        J[~np.isfinite(J).all(axis=(1, 2))] = 0.0
        A = np.einsum('bik,bil->bkl', J, J); g = np.einsum('bik,bi->bk', J, Ra); dg = np.einsum('bkk->bk', A)
        for trial in range(6):
            M = A + lam[ia][:, None, None] * (np.eye(n)[None] * (dg[:, :, None] + 1e-12 * (1 + dg.max(axis=1))[:, None, None]))
            try:
                dx = -np.linalg.solve(M, g[:, :, None])[:, :, 0]
            except np.linalg.LinAlgError:
                dx = np.zeros((len(ia), n))
                for b in range(len(ia)):
                    try:
                        dx[b] = -np.linalg.lstsq(M[b], g[b], rcond=None)[0]
                    except Exception:
                        pass
            xn = xa.copy(); xn[:, :n] += dx
            with np.errstate(all='ignore'):
                Rn = fun(xn); cn = np.einsum('bi,bi->b', Rn, Rn)
            cn[~np.isfinite(cn)] = np.inf
            acc = cn < cost[ia]
            sel = ia[acc]
            x[sel] = xn[acc]; R[sel] = Rn[acc]; cost[sel] = cn[acc]
            lam[sel] = np.maximum(lam[sel] / 3, 1e-15)
            rej = ~acc
            lam[ia[rej]] *= 4
            if not rej.any():
                break
            ia = ia[rej]; xa = xa[rej]; Ra = Ra[rej]; A = A[rej]; g = g[rej]; dg = dg[rej]
        active = (cost > 1e-32) & (lam < 1e12)
    return x[:, :n], cost


def main():
    av = sys.argv
    k, lin, field, sides = int(av[1]), int(av[2]), av[3], av[4]
    ns, seed = int(av[5]), int(av[6])
    sigma = float(av[7]) if len(av) > 7 else 1e-3
    B = int(av[8]) if len(av) > 8 else 500
    it1 = int(av[9]) if len(av) > 9 else 400
    it2 = int(av[10]) if len(av) > 10 else 200
    planted = int(av[11]) if len(av) > 11 else 0
    m = k - lin
    tag = (f"k{k}_l{lin}_{field}_{'m' if sides == '-' else 'mp'}_s{seed}_sg{sigma:g}" + ("_PL" if planted else ""))
    rng = np.random.default_rng(seed)
    t0 = time.time()
    X0 = starts(rng, ns, m, field)
    tau = None; Xstar = None
    if planted:
        Xstar = starts(rng, ns, m, field)
        a_s, c_s = ac_of(Xstar, field)
        tau = {s: n2.remainders(a_s, c_s, lin, s) for s in sides}
    flog, frem = build(k, lin, field, sides, sigma, tau)
    X1 = np.empty_like(X0); X2 = np.empty_like(X0)
    for b0 in range(0, ns, B):
        rows = np.arange(b0, min(ns, b0 + B))
        x1, _ = lm_rows(flog, X0[rows], rows, it1)
        x2, _ = lm_rows(frem, x1, rows, it2)
        X1[rows] = x1; X2[rows] = x2
        print(f"[{tag}] {rows[-1] + 1}/{ns} starts, {time.time() - t0:.0f}s", flush=True)
    allrows = np.arange(ns)
    res = {}
    for name, X in (('stage1', X1), ('polish', X2)):
        a, c = ac_of(X, field)
        Rr = frem(X, allrows)
        with np.errstate(all='ignore'):
            Rl = real_rows(n2.log_resid(a, c, lin, sides, field,
                                        None if tau is None else {s: tau[s] for s in sides}), field)
            gap = n2.all_gap(a, c, sides)[1]
        rr = np.abs(Rr).max(axis=1); rl = np.abs(Rl).max(axis=1)
        sep = nc.pair_sep(a, c, lin)
        mg = nc.margins(a, c) if field != 'cx' else np.full(ns, np.nan)
        imag = (np.abs(a.imag).max(axis=1) + np.abs(c.imag).max(axis=1)) if field == 'cx' else np.zeros(ns)
        extra = {}
        if planted:
            a_s, c_s = ac_of(Xstar, field)
            # distance to x* up to permutation (greedy)
            dd = np.zeros(ns)
            for i in range(ns):
                used = set()
                for j in range(m):
                    dist = np.abs(a[i] - a_s[i, j]) + np.abs(c[i] - c_s[i, j])
                    for l in used:
                        dist[l] = np.inf
                    l = int(np.argmin(dist)); used.add(l); dd[i] = max(dd[i], dist[l])
            extra['dstar'] = dd
        np.savez_compressed(f"ns_out/{tag}_{name}.npz", x=X, rr=rr, rl=rl, sep=sep, gap=gap, mg=mg, imag=imag,
                            a=a, c=c, **extra)
        q = np.quantile(rr, [0, 0.01, 0.1, 0.5, 0.9])
        print(f"== {tag} [{name}] remainder residual quantiles (0,1,10,50,90%): " + " ".join(f"{z:.2e}" for z in q))
        q = np.quantile(rl, [0, 0.01, 0.1, 0.5, 0.9])
        print(f"   log residual quantiles: " + " ".join(f"{z:.2e}" for z in q))
        good = (sep >= 1e-3) & (gap >= 1e-3)
        print(f"   non-degenerate (sep >= 1e-3 and all special-point gaps >= 1e-3): {good.sum()} of {ns}")
        for thr in (1e-12, 1e-9, 1e-6, 1e-4, 1e-2):
            z = rr < thr
            print(f"   rem residual < {thr:g}: {z.sum()}; non-degenerate: {(z & good).sum()}"
                  + (f"; equal to x*: {(z & good & (extra['dstar'] < 1e-6)).sum()}" if planted else ""))
        for dl in (1e-3, 1e-2, 0.05, 0.1, 0.3, 1.0):
            ok = (sep >= dl) & (gap >= dl)
            if ok.any():
                i = np.argmin(np.where(ok, rr, np.inf))
                print(f"   sep, gap >= {dl:<5g}: n = {ok.sum():6d}, min rem residual {rr[i]:.3e} (log {rl[i]:.2e}, sep "
                      f"{sep[i]:.2e}, gap {gap[i]:.2e}, margin {mg[i]:.2e}, imag {imag[i]:.2e})")
        order = np.argsort(np.where(good, rr, np.inf))[:3]
        for i in order:
            if good[i]:
                print(f"   best: rem {rr[i]:.3e} log {rl[i]:.3e} sep {sep[i]:.3e} gap {gap[i]:.3e} margin {mg[i]:.3e}"
                      f"\n      a = {np.array2string(a[i], precision=8, max_line_width=400)}"
                      f"\n      c = {np.array2string(c[i], precision=8, max_line_width=400)}")
    print(f"total time {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
