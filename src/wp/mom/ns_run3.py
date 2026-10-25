#!/usr/bin/env python3
"""ns line, driver v3 (polynomial-identity + variable projection, ns_core3).  Usage:
  ns_run3.py k lin field sides nstarts seed [batch=500] [iters=500] [planted=0] [sigma=1e-3]
field: real | cx | box ; sides: -+ | - .  planted = 1: each start gets its own random x* (same field), the system
is shifted by the fixed polynomial defect of x*, so x* is an exact zero (recall control).
sigma > 0 (default 1e-3): deflation 1 + sigma sum_i 1/|a_i|^2 against the degenerate attractors of this
formulation: the cubic coefficient a_i -> 0 (member degenerates to a quadratic, or to the linear 2t + 1 when also
c_i -> 0; the identity stays solvable in the limit, e.g. (4t^3 + 4t^2 + 2t + 1) + 4 x (2t + 1) solves the k = 5
"-" half because (2t+1)^4 == 1 mod 2t^2 + 2t + 1).  In the real domain c^2 < 8a, so only the linear collapse occurs.
Reports: identity residual (projected, weighted, max-abs), the independent per-member remainder residual (ns_core),
separation of members, gaps between special points, reality margin.
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ns_core as nc
import ns_core2 as n2
import ns_core3 as n3
from ns_run2 import _lm_fixed_last

AMIN = float(os.environ.get('NS_AMIN', '0.01'))   # lower bound for the cubic coefficient a_i (real / box fields)
DPOW = float(os.environ.get('NS_DPOW', '2'))      # deflation power for |a_i| (complex field)
WREM = float(os.environ.get('NS_WREM', '0'))      # > 0: append WREM x (per-member remainder residual, ns_core) to the
                                                   # LM residual (kills the identity residual's zeros at infinity)
RINF = float(os.environ.get('NS_RINF', '0'))      # > 0: deflation factor 1 + sum(|a|^2 + |c|^2)/RINF^2 vs escape to infinity
RMAX = np.sqrt(8 / AMIN) - 0.5


def sig(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -40, 40)))


WMODE = os.environ.get('NS_WMODE', 'free')        # 'free': W any (fixed lead); 'c': Laurent shape + universal constraints
RPOW = float(os.environ.get('NS_RPOW', '1'))
PAIR = float(os.environ.get('NS_PAIR', '0'))      # > 0: deflation 1 + PAIR sum_{i<j} 1/d_ij^2 (members in (a, c))
BDEF = float(os.environ.get('NS_BDEF', '0'))      # > 0: deflation 1 + BDEF sum_i 1/margin_i^2 against the reality boundary
                                                   # (real fields; margin_i = (sqrt(8a) - max(|c|, |a-c|))/sqrt(8a))      # power of the escape-to-infinity deflation factor
PLUSJ = int(os.environ.get('NS_PLUSJ', '-1'))     # >= 0: partial system: '-' for all members, '+' for the first PLUSJ
SYMQ = int(os.environ.get('NS_SYMQ', '1'))        # field 'sym': number of self-reflected members (u = 1/2)


def ac_of(x, field):
    if field == 'sym':
        # reflection-symmetric family: p pairs (u, v), (1 - u, v) and q = SYMQ members with u = 1/2
        q = SYMQ; p = (x.shape[1] - q) // 2
        rp = RMAX * sig(x[:, :p]); vp = (rp + 0.5) ** 2 / 8; up = 0.5 + rp * np.tanh(x[:, p:2 * p])
        rq = RMAX * sig(x[:, 2 * p:]); vq = (rq + 0.5) ** 2 / 8; uq = 0.5 + 0 * rq
        u = np.concatenate([up, 1 - up, uq], axis=1); v = np.concatenate([vp, vp, vq], axis=1)
        return 1.0 / v + 0j, u / v + 0j
    return ac_of0(x, field)


def ac_of0(x, field):
    """real: exact reality domain (|u - 1/2| + 1/2)^2 < 8v intersected with a > AMIN (v < 1/AMIN):
       r = RMAX sig(rho), v = (r + 1/2)^2/8, u = 1/2 + r tanh(z), a = 1/v, c = u/v.
       box: a = AMIN + (40 - AMIN) sig(alpha), c free.   cx: complex (a, c) free."""
    m = x.shape[1] // 2
    if field == 'real':
        r = RMAX * sig(x[:, :m]); v = (r + 0.5) ** 2 / 8; u = 0.5 + r * np.tanh(x[:, m:])
        return 1.0 / v + 0j, u / v + 0j
    if field == 'box':
        return AMIN + (40.0 - AMIN) * sig(x[:, :m]) + 0j, x[:, m:] + 0j
    return nc.cx_unpack(x)


def logit(p):
    return np.log(p / (1 - p))


def starts(rng, B, m, field):
    if field == 'sym':
        q = SYMQ; p = (m - q) // 2
        a = np.exp(rng.uniform(np.log(max(0.05, 2 * AMIN)), np.log(31.5), (B, p + q)))
        r = np.sqrt(8 / a) - 0.5
        return np.concatenate([logit(r[:, :p] / RMAX), rng.normal(0, 1.0, (B, p)), logit(r[:, p:] / RMAX)], axis=1)
    a = np.exp(rng.uniform(np.log(max(0.05, 2 * AMIN)), np.log(31.5), (B, m)))
    if field == 'real':
        r = np.sqrt(8 / a) - 0.5
        return np.concatenate([logit(r / RMAX), rng.normal(0, 1.0, (B, m))], axis=1)
    c = rng.uniform(-1, 1, (B, m)) * np.sqrt(8 * a) * 0.9 + 0.5 * a
    if field == 'box':
        return np.concatenate([logit((a - AMIN) / (40 - AMIN)), c], axis=1)
    sa = np.exp(rng.uniform(np.log(0.01), np.log(2), (B, 1)))
    ai = rng.normal(0, 1, (B, m)) * sa * a; ci = rng.normal(0, 1, (B, m)) * sa * (np.abs(c) + 1)
    return np.concatenate([a, ai, c, ci], axis=1)


def main():
    av = sys.argv
    k, lin, field, sides = int(av[1]), int(av[2]), av[3], av[4]
    ns, seed = int(av[5]), int(av[6])
    B = int(av[7]) if len(av) > 7 else 500
    iters = int(av[8]) if len(av) > 8 else 500
    planted = int(av[9]) if len(av) > 9 else 0
    sigma = float(av[10]) if len(av) > 10 else 1e-3
    m = k - lin
    real = field != 'cx'
    print(f"AMIN = {AMIN}, DPOW = {DPOW}, RINF = {RINF}^{RPOW}, WREM = {WREM}, WMODE = {WMODE}, BDEF = {BDEF}, PAIR = {PAIR}" + (f", SYMQ = {SYMQ}" if field == 'sym' else ""))
    if field == 'sym' and (m - SYMQ) % 2:
        raise SystemExit("m - SYMQ must be even")
    tag = f"v3_k{k}_l{lin}_{field}{SYMQ if field == 'sym' else ''}_" + (f"J{PLUSJ}_" if PLUSJ >= 0 else "") + f"{'m' if sides == '-' else 'mp'}_s{seed}" + ("_PL" if planted else "") + \
          f"_sg{sigma:g}" + (f"_w{WREM:g}" if WREM > 0 else "") + ("_Wc" if WMODE == 'c' else "") + (f"_b{BDEF:g}" if BDEF > 0 else "") + (f"_p{PAIR:g}" if PAIR > 0 else "")
    y = n3.eval_points(k)
    rng = np.random.default_rng(seed)
    t0 = time.time()
    X0 = starts(rng, ns, m, field)
    defect = None; Xs = None; tau = None
    if planted:
        Xs = starts(rng, ns, m, field)
        a_s, c_s = ac_of(Xs, field)
        defect = n3.defects_c(a_s, c_s, lin, y, sides, real) if WMODE == 'c' else n3.defects(a_s, c_s, lin, y, sides, real)
        tau = {s: n2.remainders(a_s, c_s, lin, s) for s in sides}

    def remrows(a, c, rows):
        Rm = n2.rem_resid(a, c, lin, sides, None if tau is None else {s: (tau[s][0][rows], tau[s][1][rows]) for s in sides})
        if PLUSJ >= 0:   # columns: '-' r0 (m), '-' r1 (m), '+' r0 (m), '+' r1 (m); keep '+' for the first PLUSJ only
            mm = a.shape[1]
            Rm = np.concatenate([Rm[:, :2 * mm], Rm[:, 2 * mm:2 * mm + PLUSJ], Rm[:, 3 * mm:3 * mm + PLUSJ]], axis=1)
        Rm = np.concatenate([Rm.real, Rm.imag], axis=1) if field == 'cx' else Rm.real
        return np.where(np.isfinite(Rm), Rm, 1e3)

    def fun(z):
        rows = z[:, -1].astype(int); x = z[:, :-1]
        a, c = ac_of(x, field)
        with np.errstate(all='ignore'):
            if PLUSJ >= 0:
                R = n3.resid_partial(a, c, lin, y, m, PLUSJ, real)
            elif WMODE == 'c':
                R = n3.resid_c_def(a, c, lin, y, sides, real, None if defect is None else {s: defect[s][rows] for s in sides})
            else:
                R = n3.resid(a, c, lin, y, sides, real, None if defect is None else {s: defect[s][rows] for s in sides})
            if WREM > 0:
                R = np.concatenate([R, WREM * remrows(a, c, rows)], axis=1)
            if sigma > 0 or BDEF > 0 or PAIR > 0:
                D = 1.0 + sigma * (1.0 / (np.abs(a) ** DPOW + 1e-300)).sum(axis=1)
                if BDEF > 0 and field in ('real', 'sym'):
                    ar = np.maximum(a.real, 1e-300); s8 = np.sqrt(8 * ar)
                    mgi = (s8 - np.maximum(np.abs(c.real), np.abs(ar - c.real))) / s8
                    D = D * (1.0 + BDEF * (1.0 / np.maximum(mgi, 1e-150) ** 2).sum(axis=1))
                if PAIR > 0:
                    D = D * (1.0 + PAIR * n2.pair_terms(a, c, lin))
                if RINF > 0:
                    D = D * (1.0 + ((np.abs(a) ** 2 + np.abs(c) ** 2).sum(axis=1)) / RINF ** 2) ** RPOW
                R = R * np.where(np.isfinite(D), D, 1e300)[:, None]
        return R

    X = np.empty_like(X0); cost = np.empty(ns)
    for b0 in range(0, ns, B):
        rows = np.arange(b0, min(ns, b0 + B))
        z = np.concatenate([X0[rows], rows[:, None].astype(float)], axis=1)
        x, cst = _lm_fixed_last(fun, z, iters)
        X[rows] = x; cost[rows] = cst
        print(f"[{tag}] {rows[-1] + 1}/{ns}, {time.time() - t0:.0f}s", flush=True)
    a, c = ac_of(X, field)
    with np.errstate(all='ignore'):
        if PLUSJ >= 0:
            R3 = n3.resid_partial(a, c, lin, y, m, PLUSJ, real)
        elif WMODE == 'c':
            R3 = n3.resid_c_def(a, c, lin, y, sides, real, defect)
            R3f = n3.resid(a, c, lin, y, sides, real, None) if defect is None else None
            if R3f is not None:
                rrf = np.abs(R3f).max(axis=1)
        else:
            R3 = n3.resid(a, c, lin, y, sides, real, defect)
        rr = np.abs(R3).max(axis=1)
        rrc = rr.copy()
        if WMODE == 'c' and defect is None and PLUSJ < 0:
            rr = rrf     # report the (better normalised) free-W identity residual; same zero set; rrc saved too
        rem = np.abs(remrows(a, c, np.arange(ns))).max(axis=1)
        if planted:
            a_s, c_s = ac_of(Xs, field)
            dstar = np.zeros(ns)
            for i in range(ns):
                used = set()
                for j in range(m):
                    d = np.abs(a[i] - a_s[i, j]) + np.abs(c[i] - c_s[i, j])
                    for l in used:
                        d[l] = np.inf
                    l = int(np.argmin(d)); used.add(l); dstar[i] = max(dstar[i], d[l])
        else:
            dstar = np.full(ns, np.nan)
        gap = n2.all_gap(a, c, sides)[1]
    sep = nc.pair_sep(a, c, lin)
    mg = nc.margins(a, c) if real else np.full(ns, np.nan)
    mina = np.abs(a).min(axis=1)
    imag = (np.abs(a.imag).max(axis=1) + np.abs(c.imag).max(axis=1)) if field == 'cx' else np.zeros(ns)
    os.makedirs('ns_out', exist_ok=True)
    np.savez_compressed(f"ns_out/{tag}.npz", x=X, rr=rr, rrc=rrc, rem=rem, sep=sep, gap=gap, mg=mg, mina=mina, imag=imag,
                        dstar=dstar, a=a, c=c, X0=X0, Xs=Xs if Xs is not None else np.zeros(0))
    q = np.quantile(rr, [0, 0.01, 0.1, 0.5, 0.9])
    print(f"== {tag}: identity residual quantiles (0,1,10,50,90%): " + " ".join(f"{z:.2e}" for z in q))
    if WMODE == 'c':
        q = np.quantile(rrc, [0, 0.01, 0.1, 0.5, 0.9])
        print(f"   constrained-W identity residual quantiles: " + " ".join(f"{z:.2e}" for z in q))
    q = np.nanquantile(rem, [0, 0.01, 0.1, 0.5, 0.9])
    print(f"   remainder residual quantiles: " + " ".join(f"{z:.2e}" for z in q))
    both = np.maximum(rr, rem)
    q = np.nanquantile(both, [0, 0.01, 0.1, 0.5, 0.9])
    print(f"   max(identity, remainder) quantiles: " + " ".join(f"{z:.2e}" for z in q))
    good = (sep >= 1e-3) & (gap >= 1e-3) & (mina >= 1e-3)
    print(f"   non-degenerate (sep, special-point gaps, min|a| >= 1e-3): {good.sum()} of {ns}")
    for thr in (1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2):
        z = rr < thr
        s = (f"   identity residual < {thr:g}: {z.sum()}; non-degenerate: {(z & good).sum()}; and remainder < {thr:g}: "
             f"{(z & good & (rem < thr)).sum()}")
        if planted:
            s += f"; = x* (dist < 1e-6): {(z & (dstar < 1e-6)).sum()}"
        print(s)
    for dl in (1e-3, 1e-2, 0.05, 0.1, 0.3):
        ok = (sep >= dl) & (gap >= dl) & (mina >= dl)
        if ok.any():
            i = np.argmin(np.where(ok, both, np.inf))
            print(f"   sep, gap, min|a| >= {dl:<5g}: n = {ok.sum():6d}, min max(id, rem) {both[i]:.3e}: identity {rr[i]:.3e} (rem {rem[i]:.2e},"
                  f" sep {sep[i]:.2e}, gap {gap[i]:.2e}, min|a| {mina[i]:.2e}, margin {mg[i]:.2e}, imag {imag[i]:.2e})")
    order = np.argsort(np.where(good, both, np.inf))[:4]
    for i in order:
        if good[i]:
            print(f"   best: id {rr[i]:.3e} rem {rem[i]:.3e} sep {sep[i]:.3e} gap {gap[i]:.3e} margin {mg[i]:.3e}"
                  f"\n      a = {np.array2string(a[i], precision=10, max_line_width=500)}"
                  f"\n      c = {np.array2string(c[i], precision=10, max_line_width=500)}")
    print(f"total time {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
