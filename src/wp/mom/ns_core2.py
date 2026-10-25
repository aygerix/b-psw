#!/usr/bin/env python3
"""ns line, core v2: LOG residual for the global phase (gradient does not saturate when |P_i| >> 1) and the exact
REMAINDER residual (ns_core) for polishing / verification.

Condition (member i, side s): P_i == tau0 + tau1 t  mod  m_{i,s}  (m_{i,-} = g_i/a_i, m_{i,+} = h_i/a_i);
the Williams family is tau = (1, 0) everywhere; a PLANTED control uses tau = remainder of P_i(x*) at a random x*.
Log residual at a root w of m_{i,s}: Q = P_i(w)/(tau0 + tau1 w); components log|Q|, Re(Q/|Q|) - 1, Im(Q/|Q|)
(continuous, zero only at Q = 1).  Real field: only the upper root (conjugate gives the conjugate equation).
Complex field: both roots.  Evaluation at roots misses the derivative condition at a double root, so the deflation
repels ALL coincidences among special points (same member included) and the shared points 0, -1.
"""
import numpy as np
import ns_core as nc


def fval(a, c, t):
    """f_j(t) for t of shape (B, n) and members (B, m): returns (B, n, m)"""
    T = t[:, :, None]; A = a[:, None, :]; C = c[:, None, :]
    return 2 * T + 1 + T * (T + 1) * (A * T + C)


def mod_pq(a, c, side):
    if side == '-':
        return (a + c) / a, (c + 2) / a
    return c / a, 2 / a


def roots_pq(p, q):
    d = np.sqrt(p * p - 4 * q + 0j)
    r1 = (-p + d) / 2; r2 = (-p - d) / 2
    # order: upper root first (for real data Im r1 >= 0 picks the upper one)
    sw = r1.imag < r2.imag
    return np.where(sw, r2, r1), np.where(sw, r1, r2)


def log_resid(a, c, lin, sides, field, tau=None):
    """tau: dict side -> (tau0, tau1) arrays (B, m) or None for the Williams targets (1, 0)."""
    B, m = a.shape
    out = []
    idx = np.arange(m)
    for side in sides:
        p, q = mod_pq(a, c, side)
        w1, w2 = roots_pq(p, q)
        roots = [w1] if field != 'cx' else [w1, w2]
        for w in roots:
            fv = fval(a, c, w)                                  # (B, i, j): f_j at root of member i
            fv[:, idx, idx] = 1.0
            lmag = np.log(np.abs(fv)).sum(axis=2)
            ph = (fv / np.abs(fv)).prod(axis=2)
            if lin:
                l = 2 * w + 1
                lmag = lmag + np.log(np.abs(l)); ph = ph * (l / np.abs(l))
            if tau is not None:
                T = tau[side][0] + tau[side][1] * w
                lmag = lmag - np.log(np.abs(T)); ph = ph * (np.abs(T) / T)
            out += [lmag, ph.real - 1.0, ph.imag]
    return np.concatenate(out, axis=1)


def rem_resid(a, c, lin, sides, tau=None):
    """normalised remainder residual (ns_core.side_resid) with optional targets"""
    out = []
    for side in sides:
        e0, e1, p, q = nc.side_resid(a, c, lin, side, raw=True)
        if tau is not None:
            e0 = e0 + 1.0 - tau[side][0]; e1 = e1 - tau[side][1]
        # scale as in ns_core
        x0, x1 = nc.reduce_members(a, c, p, q)
        m = a.shape[1]; idx = np.arange(m)
        rho = np.sqrt(np.abs(q))
        lt = np.log1p(np.abs(x0) + np.abs(x1) * rho[:, :, None]); lt[:, idx, idx] = 0.0
        ls = lt.sum(axis=2)
        if lin:
            ls = ls + np.log1p(1.0 + 2.0 * rho)
        S = np.exp(ls)
        out += [e0 / S, e1 * rho / S]
    return np.concatenate(out, axis=1)


def remainders(a, c, lin, side):
    """raw remainder (r0, r1) of P_i mod m_{i,side} (so the target is (1, 0))"""
    e0, e1, p, q = nc.side_resid(a, c, lin, side, raw=True)
    return e0 + 1.0, e1


def all_gap(a, c, sides):
    """sum 1/|w - w'|^2 over all pairs of special points (both roots of every member, plus 0 and -1; the pair
    (0, -1) excluded), per side; and the min gap."""
    B, m = a.shape
    n = 2 * m + 2
    iu = np.triu_indices(n, 1)
    keep = ~((iu[0] == 2 * m) & (iu[1] == 2 * m + 1))
    tot = np.zeros(B); gmin = np.full(B, np.inf)
    for side in sides:
        p, q = mod_pq(a, c, side)
        w1, w2 = roots_pq(p, q)
        w = np.concatenate([w1, w2, np.zeros((B, 1)), -np.ones((B, 1))], axis=1)
        d2 = np.abs(w[:, :, None] - w[:, None, :]) ** 2
        d2 = d2[:, iu[0], iu[1]][:, keep]
        tot += (1.0 / (d2 + 1e-300)).sum(axis=1)
        gmin = np.minimum(gmin, np.sqrt(d2.min(axis=1)))
    return tot, gmin


def pair_terms(a, c, lin):
    B, m = a.shape
    d = np.abs(a[:, :, None] - a[:, None, :]) ** 2 + np.abs(c[:, :, None] - c[:, None, :]) ** 2
    iu = np.triu_indices(m, 1)
    s = (1.0 / (d[:, iu[0], iu[1]] + 1e-300)).sum(axis=1)
    if lin:
        s = s + (1.0 / (np.abs(a) ** 2 + np.abs(c) ** 2 + 1e-300)).sum(axis=1)
    return s


def deflation(a, c, lin, sigma, sides):
    return 1.0 + sigma * (pair_terms(a, c, lin) + all_gap(a, c, sides)[0])
