#!/usr/bin/env python3
"""ns line (direct numerical search for real cubic Williams families), core routines.

Normal form (Lemma 24c): f_i = 2t + 1 + t(t+1)(a_i t + c_i) = a t^3 + (a+c) t^2 + (c+2) t + 1,
g_i = (f_i - 1)/t = a t^2 + (a+c) t + c + 2,   h_i = (f_i + 1)/(t+1) = a t^2 + c t + 2,
optional linear member 2t + 1 (its own two conditions are automatic: F(0) = 1, F(-1) = -1 for k odd).
Family condition (k odd): P_i := prod_{j != i} f_j == 1 mod g_i and == 1 mod h_i, for every cubic i.
(F - 1 = f_i P_i - 1 == P_i - 1 mod (f_i - 1); F + 1 == 1 - P_i mod (f_i + 1); t, t+1 parts automatic.)

Arithmetic in C[t]/(t^2 + p t + q): elements x0 + x1 t, (x0 + x1 t)(y0 + y1 t) = (x0 y0 - q x1 y1) + (x0 y1 + x1 y0 - p x1 y1) t.
"-" modulus g_i/a_i: p = (a+c)/a = 1 + u, q = (c+2)/a = u + 2v.   "+" modulus h_i/a_i: p = c/a = u, q = 2/a = 2v.
f_j mod (t^2 + p t + q): t^3 = (p^2 - q) t + p q, so
    x1 = a_j (p^2 - q) - (a_j + c_j) p + c_j + 2,   x0 = a_j p q - (a_j + c_j) q + 1.
Normalised residual per (member, side): r0 = (P0 - 1)/S, r1 = P1 rho/S, rho = sqrt|q| (~ |root|),
S = prod_{j != i} (1 + |x0_j| + |x1_j| rho)  (scale of |P_i(root)|; residual is scale free, <~ 1).
Everything is batched: a, c have shape (B, m), complex128.
"""
import numpy as np


def reduce_members(a, c, p, q):
    """x0, x1 of shape (B, m_i, m_j): f_j reduced mod t^2 + p_i t + q_i"""
    P = p[:, :, None]; Q = q[:, :, None]; A = a[:, None, :]; C = c[:, None, :]
    x1 = A * (P * P - Q) - (A + C) * P + C + 2
    x0 = A * P * Q - (A + C) * Q + 1
    return x0, x1


def side_resid(a, c, lin, side, raw=False):
    """normalised remainders (r0, r1), each (B, m) complex, for side '-' (mod g_i) or '+' (mod h_i)"""
    if side == '-':
        p = (a + c) / a; q = (c + 2) / a
    else:
        p = c / a; q = 2 / a
    B, m = a.shape
    x0, x1 = reduce_members(a, c, p, q)
    idx = np.arange(m)
    x0[:, idx, idx] = 1.0; x1[:, idx, idx] = 0.0
    rho = np.sqrt(np.abs(q))
    lt = np.log1p(np.abs(x0) + np.abs(x1) * rho[:, :, None])
    lt[:, idx, idx] = 0.0
    lscale = lt.sum(axis=2)
    P0 = np.ones((B, m), dtype=complex); P1 = np.zeros((B, m), dtype=complex)
    for j in range(m):
        y0 = x0[:, :, j]; y1 = x1[:, :, j]
        P0, P1 = P0 * y0 - q * P1 * y1, P0 * y1 + P1 * y0 - p * P1 * y1
    if lin:
        P0, P1 = P0 * 1.0 - q * P1 * 2.0, P0 * 2.0 + P1 * 1.0 - p * P1 * 2.0
        lscale = lscale + np.log1p(1.0 + 2.0 * rho)
    if raw:
        return P0 - 1.0, P1, p, q
    S = np.exp(lscale)
    return (P0 - 1.0) / S, P1 * rho / S


def full_resid(a, c, lin, sides='-+'):
    out = []
    for sd in sides:
        r0, r1 = side_resid(a, c, lin, sd)
        out += [r0, r1]
    return np.concatenate(out, axis=1)  # (B, 2 m len(sides)) complex


# ---------------- parametrisations ----------------
def real_unpack(x):
    """x (B, 2m) real -> (u, v); exact parametrisation of the reality domain (|u-1/2|+1/2)^2 < 8v (m9 line)"""
    m = x.shape[1] // 2
    r = np.exp(np.clip(x[:, :m], -40, 40))
    v = (r + 0.5) ** 2 / 8
    u = 0.5 + r * np.tanh(x[:, m:])
    return u, v


def real_pack(u, v):
    r = np.sqrt(8 * v) - 0.5
    z = np.arctanh(np.clip((u - 0.5) / r, -1 + 1e-15, 1 - 1e-15))
    return np.concatenate([np.log(r), z], axis=1)


def uv_to_ac(u, v):
    return 1.0 / v, u / v


def ac_to_uv(a, c):
    return c / a, 1.0 / a


def free_unpack(x):
    """x (B, 2m) real -> (a, c) real, no domain restriction"""
    m = x.shape[1] // 2
    return x[:, :m].astype(complex), x[:, m:].astype(complex)


def cx_unpack(x):
    """x (B, 4m) real -> complex (a, c)"""
    m = x.shape[1] // 4
    return x[:, :m] + 1j * x[:, m:2 * m], x[:, 2 * m:3 * m] + 1j * x[:, 3 * m:]


# ---------------- diagnostics ----------------
def pair_sep(a, c, lin):
    """min pairwise distance in (a, c) coefficient space (incl. distance to the linear member (0, 0) if lin)"""
    B, m = a.shape
    d = np.abs(a[:, :, None] - a[:, None, :]) ** 2 + np.abs(c[:, :, None] - c[:, None, :]) ** 2
    idx = np.arange(m)
    d[:, idx, idx] = np.inf
    dm = np.sqrt(d.min(axis=(1, 2)))
    if lin:
        dm = np.minimum(dm, np.sqrt((np.abs(a) ** 2 + np.abs(c) ** 2).min(axis=1)))
    return dm


def special_roots(a, c, side):
    """roots of g_i (side '-') or h_i (side '+'), shape (B, 2m) complex, owner index of each root"""
    if side == '-':
        p = (a + c) / a; q = (c + 2) / a
    else:
        p = c / a; q = 2 / a
    dsc = np.sqrt(p * p - 4 * q + 0j)
    return np.concatenate([(-p + dsc) / 2, (-p - dsc) / 2], axis=1)


def root_gap_terms(a, c, sides='-+'):
    """sum over pairs of special points with DIFFERENT owners (members, and the shared points 0, -1) of 1/|w - w'|^2,
    separately for the '-' and '+' points; also the minimum such gap.  Coincidences of this kind are spurious zeros
    of the remainder residual (they hide a derivative condition, or put a root where P_i = 1 automatically)."""
    B, m = a.shape
    own = np.concatenate([np.arange(m), np.arange(m), [m, m + 1]])
    diff = own[:, None] != own[None, :]
    diff[2 * m, 2 * m + 1] = diff[2 * m + 1, 2 * m] = False
    iu = np.triu_indices(2 * m + 2, 1)
    keep = diff[iu]
    tot = np.zeros(B); gmin = np.full(B, np.inf)
    for side in sides:
        w = special_roots(a, c, side)
        w = np.concatenate([w, np.zeros((B, 1)), -np.ones((B, 1))], axis=1)
        d2 = np.abs(w[:, :, None] - w[:, None, :]) ** 2
        d2 = d2[:, iu[0], iu[1]][:, keep]
        tot += (1.0 / (d2 + 1e-300)).sum(axis=1)
        gmin = np.minimum(gmin, np.sqrt(d2.min(axis=1)))
    return tot, gmin


def deflation(a, c, lin, sigma, roots=True, sides='-+'):
    B, m = a.shape
    d = np.abs(a[:, :, None] - a[:, None, :]) ** 2 + np.abs(c[:, :, None] - c[:, None, :]) ** 2
    iu = np.triu_indices(m, 1)
    D = 1.0 + sigma * (1.0 / (d[:, iu[0], iu[1]] + 1e-300)).sum(axis=1)
    if lin:
        D = D + sigma * (1.0 / (np.abs(a) ** 2 + np.abs(c) ** 2 + 1e-300)).sum(axis=1)
    if roots:
        D = D + sigma * root_gap_terms(a, c, sides)[0]
    return D


def margins(a, c):
    """real domain margins: min_i (sqrt(8a) - max(|c|, |a - c|)) / sqrt(8a)  (> 0 iff g_i, h_i have no real root)"""
    a = a.real; c = c.real
    s8 = np.sqrt(8 * np.maximum(a, 1e-300))
    return ((s8 - np.maximum(np.abs(c), np.abs(a - c))) / s8).min(axis=1)


# ---------------- batched Levenberg-Marquardt ----------------
def batched_lm(fun, x0, iters=400, lam0=1e-2, h=1e-7, tol=1e-30, verbose=False):
    """fun: (B, n) real -> (B, N) real. Minimises ||fun||^2 per row independently. Returns x, cost history."""
    x = x0.copy()
    B, n = x.shape
    R = fun(x)
    cost = np.einsum('bi,bi->b', R, R)
    lam = np.full(B, lam0)
    active = np.isfinite(cost)
    cost[~active] = np.inf
    for it in range(iters):
        ia = np.nonzero(active)[0]
        if len(ia) == 0:
            break
        xa = x[ia]; Ra = R[ia]
        J = np.empty((len(ia), Ra.shape[1], n))
        for col in range(n):
            hh = h * (1 + np.abs(xa[:, col]))
            x2 = xa.copy(); x2[:, col] += hh
            J[:, :, col] = (fun(x2) - Ra) / hh[:, None]
        bad = ~np.isfinite(J).all(axis=(1, 2))
        J[bad] = 0.0
        A = np.einsum('bik,bil->bkl', J, J)
        g = np.einsum('bik,bi->bk', J, Ra)
        dg = np.einsum('bkk->bk', A)
        for trial in range(6):
            M = A + (lam[ia][:, None, None]) * (np.eye(n)[None] * (dg[:, :, None] + 1e-12 * (1 + dg.max(axis=1))[:, None, None]))
            try:
                dx = -np.linalg.solve(M, g[:, :, None])[:, :, 0]
            except np.linalg.LinAlgError:
                dx = np.zeros_like(xa)
                for b in range(len(ia)):
                    try:
                        dx[b] = -np.linalg.lstsq(M[b], g[b], rcond=None)[0]
                    except Exception:
                        pass
            xn = xa + dx
            with np.errstate(all='ignore'):
                Rn = fun(xn)
                cn = np.einsum('bi,bi->b', Rn, Rn)
            cn[~np.isfinite(cn)] = np.inf
            acc = cn < cost[ia]
            # accept improved rows
            sel = ia[acc]
            x[sel] = xn[acc]; R[sel] = Rn[acc]; cost[sel] = cn[acc]
            lam[sel] = np.maximum(lam[sel] / 3, 1e-15)
            rej = ~acc
            lam[ia[rej]] *= 4
            if not rej.any():
                break
            # retry only rejected rows
            keep = rej
            ia = ia[keep]; xa = xa[keep]; Ra = Ra[keep]; A = A[keep]; g = g[keep]; dg = dg[keep]
            if len(ia) == 0:
                break
        active = (cost > tol) & (lam < 1e12)
        if verbose and it % 50 == 0:
            print(it, np.sort(cost)[:5], flush=True)
    return x, cost
