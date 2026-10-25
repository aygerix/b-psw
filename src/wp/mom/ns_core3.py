#!/usr/bin/env python3
"""ns line, core v3: POLYNOMIAL-IDENTITY formulation with variable projection.

"-" side:  F - 1 = t * G * W_-,   G = prod_i g_i (cubic members only),  W_- of degree k-1, leading coeff 2^lin.
"+" side:  F + 1 = (t+1) * H * W_+, H = prod_i h_i,                       W_+ of degree k-1, leading coeff 2^lin.
(lead: F has lead 2^lin prod a_i, G and H have lead prod a_i.)  These identities are EXACTLY the family conditions
(t G | F - 1 contains every g_i | (F-1)/t; conversely the g_i are pairwise coprime for a non-degenerate real family),
and they are STRONGER than the per-member conditions when members coincide: the diagonal f_1 = ... = f_k is NOT a
solution (g^2 would have to divide (F - 1)/t), so no deflation is needed against coincident members.
The difference D_-(y) = F - 1 - y G W_- has D_-(0) = 0 and degree <= 3k - 1 (top terms cancel), so it vanishes iff it
vanishes at >= 3k - 1 distinct nonzero points.  We evaluate at points on two circles, weight by 1/(1 + |F(y)|), and
eliminate the free coefficients of W (linear) by least squares (variable projection): the residual is a function
of the members only.  Planted control: subtract a fixed defect D*(y) computed at a random x*.
"""
import numpy as np


def eval_points(k, radii=(0.8, 1.6), extra=3):
    L = 3 * k + extra
    pts = [r * np.exp(2j * np.pi * (np.arange(L) + 0.5) / L) for r in radii]
    return np.concatenate(pts)


def member_vals(a, c, y):
    """f, g, h at points y (L,), members (B, m): each (B, L, m)"""
    Y = y[None, :, None]; A = a[:, None, :]; C = c[:, None, :]
    g = A * Y * Y + (A + C) * Y + C + 2
    h = A * Y * Y + C * Y + 2
    f = 1 + Y * g
    return f, g, h


def side_system(a, c, lin, y, side):
    """returns b (B, L) and design matrix Amat (B, L, k-1) with residual = b - Amat w, weights applied"""
    B, m = a.shape
    k = m + lin
    f, g, h = member_vals(a, c, y)
    F = f.prod(axis=2)
    if lin:
        F = F * (2 * y + 1)[None, :]
    wgt = 1.0 / (1.0 + np.abs(F))
    if side == '-':
        base = y[None, :] * g.prod(axis=2); rhs = F - 1
    else:
        base = (y + 1)[None, :] * h.prod(axis=2); rhs = F + 1
    lead = 2.0 ** lin
    b = (rhs - base * lead * y[None, :] ** (k - 1)) * wgt
    Amat = (base * wgt)[:, :, None] * (y[None, :, None] ** np.arange(k - 1)[None, None, :])
    return b, Amat, wgt


def project(b, Amat, real):
    """least-squares residual b - A w (w real if real=True), batched via QR; returns residual rows (real), w"""
    if real:
        bb = np.concatenate([b.real, b.imag], axis=1)
        AA = np.concatenate([Amat.real, Amat.imag], axis=1)
    else:
        bb = b; AA = Amat
    Q, R = np.linalg.qr(AA)
    qb = np.einsum('bli,bl->bi', Q.conj(), bb)
    res = bb - np.einsum('bli,bi->bl', Q, qb)
    try:
        w = np.linalg.solve(R, qb[:, :, None])[:, :, 0]
    except np.linalg.LinAlgError:
        w = np.full(qb.shape, np.nan)
    if not real:
        res = np.concatenate([res.real, res.imag], axis=1)
    return res, w


def resid(a, c, lin, y, sides, real, defect=None):
    out = []
    for s in sides:
        b, A, wgt = side_system(a, c, lin, y, s)
        if defect is not None:
            b = b - defect[s] * wgt
        r, _ = project(b, A, real)
        out.append(r)
    R = np.concatenate(out, axis=1)
    return np.where(np.isfinite(R), R, 1e3)


def defects(a, c, lin, y, sides, real):
    """UNweighted polynomial defect D*(y) = F - 1 - y G W* at x* (W* = LS-optimal), the planted target (B, L)"""
    d = {}
    for s in sides:
        b, A, wgt = side_system(a, c, lin, y, s)
        if real:
            bb = np.concatenate([b.real, b.imag], axis=1); AA = np.concatenate([A.real, A.imag], axis=1)
            Q, R = np.linalg.qr(AA)
            w = np.linalg.solve(R, np.einsum('bli,bl->bi', Q, bb)[:, :, None])[:, :, 0]
        else:
            Q, R = np.linalg.qr(A)
            w = np.linalg.solve(R, np.einsum('bli,bl->bi', Q.conj(), b)[:, :, None])[:, :, 0]
        d[s] = (b - np.einsum('bli,bi->bl', A, w)) / wgt
    return d


def side_system_sub(a, c, lin, y, side, nsub):
    """like side_system but the modulus is t * prod_{i < nsub} g_i (resp. (t+1) prod h_i) and W is fully free
    (degree 3k - 1 - 2 nsub, all coefficients unknown): the conditions for the first nsub members only."""
    B, m = a.shape
    k = m + lin
    f, g, h = member_vals(a, c, y)
    F = f.prod(axis=2)
    if lin:
        F = F * (2 * y + 1)[None, :]
    wgt = 1.0 / (1.0 + np.abs(F))
    if side == '-':
        base = y[None, :] * g[:, :, :nsub].prod(axis=2); rhs = F - 1
    else:
        base = (y + 1)[None, :] * h[:, :, :nsub].prod(axis=2); rhs = F + 1
    dW = 3 * m + lin - 1 - 2 * nsub
    b = rhs * wgt
    Amat = (base * wgt)[:, :, None] * (y[None, :, None] ** np.arange(dW + 1)[None, None, :])
    return b, Amat, wgt


def resid_partial(a, c, lin, y, nminus, nplus, real):
    """identity residual for: '-' conditions of the first nminus members and '+' conditions of the first nplus"""
    out = []
    for s, n in (('-', nminus), ('+', nplus)):
        if n == 0:
            continue
        b, A, wgt = side_system_sub(a, c, lin, y, s, n)
        r, _ = project(b, A, real)
        out.append(r)
    R = np.concatenate(out, axis=1)
    return np.where(np.isfinite(R), R, 1e3)


# ---------------- constrained cofactor (Laurent structure + universal constraints) ----------------
# W_- = 2^lin V,  V = t^(k-1) + (lin/2) t^(k-2) + R(t),  deg R <= k - 4,  V(-1) = 2^(1-k),  V'(-1) = (1-k) 2^(1-k)
# (V is the polynomial part of t^(k-1)[prod(1+x_i)(1+s/2)^lin - ...], x_i = O(s^3); session-8 universal constraints).
# R = R0 + sum_j beta_j (t+1)^2 t^j  (j = 0..k-6),  R0 = r0 + r1 (t+1) fixing R(-1), R'(-1).
# "+" side by reflection: F + 1 = (t+1) H W_+,  W_+(t) = 2^lin Vr(-1 - t), Vr of the same constrained shape.
def _vfix(k, lin, T):
    """fixed part of V at points T: t^(k-1) + (lin/2) t^(k-2) + R0(t)"""
    e0 = (-1.0) ** (k - 1) + (lin / 2) * (-1.0) ** (k - 2)
    e1 = (k - 1) * (-1.0) ** (k - 2) + (lin / 2) * (k - 2) * (-1.0) ** (k - 3)
    r0 = 2.0 ** (1 - k) - e0
    r1 = (1 - k) * 2.0 ** (1 - k) - e1
    return T ** (k - 1) + (lin / 2) * T ** (k - 2) + r0 + r1 * (T + 1)


def side_system_c(a, c, lin, y, side):
    B, m = a.shape
    k = m + lin
    f, g, h = member_vals(a, c, y)
    F = f.prod(axis=2)
    if lin:
        F = F * (2 * y + 1)[None, :]
    wgt = 1.0 / (1.0 + np.abs(F))
    if side == '-':
        base = y[None, :] * g.prod(axis=2); rhs = F - 1; T = y
    else:
        base = (y + 1)[None, :] * h.prod(axis=2); rhs = F + 1; T = -1 - y
    lead = 2.0 ** lin
    b = (rhs - base * lead * _vfix(k, lin, T)[None, :]) * wgt
    nb = k - 5
    if nb > 0:
        Bj = ((T + 1) ** 2)[:, None] * (T[:, None] ** np.arange(nb)[None, :])     # (L, nb)
        Amat = (base * wgt * lead)[:, :, None] * Bj[None, :, :]
    else:
        Amat = np.zeros((B, len(y), 0), dtype=complex)
    return b, Amat, wgt


def resid_c(a, c, lin, y, sides, real):
    out = []
    for s in sides:
        b, A, wgt = side_system_c(a, c, lin, y, s)
        if A.shape[2] > 0:
            r, _ = project(b, A, real)
        else:
            r = np.concatenate([b.real, b.imag], axis=1)
        out.append(r)
    R = np.concatenate(out, axis=1)
    return np.where(np.isfinite(R), R, 1e3)


def defects_c(a, c, lin, y, sides, real):
    """planted target for the constrained formulation: UNweighted defect D* = rhs - base W*_c at x*"""
    d = {}
    for s in sides:
        b, A, wgt = side_system_c(a, c, lin, y, s)
        if A.shape[2] > 0:
            if real:
                bb = np.concatenate([b.real, b.imag], axis=1); AA = np.concatenate([A.real, A.imag], axis=1)
                Q, R = np.linalg.qr(AA)
                w = np.linalg.solve(R, np.einsum('bli,bl->bi', Q, bb)[:, :, None])[:, :, 0]
            else:
                Q, R = np.linalg.qr(A)
                w = np.linalg.solve(R, np.einsum('bli,bl->bi', Q.conj(), b)[:, :, None])[:, :, 0]
            d[s] = (b - np.einsum('bli,bi->bl', A, w)) / wgt
        else:
            d[s] = b / wgt
    return d


def resid_c_def(a, c, lin, y, sides, real, defect):
    out = []
    for s in sides:
        b, A, wgt = side_system_c(a, c, lin, y, s)
        if defect is not None:
            b = b - defect[s] * wgt
        if A.shape[2] > 0:
            r, _ = project(b, A, real)
        else:
            r = np.concatenate([b.real, b.imag], axis=1) if real else np.concatenate([b.real, b.imag], axis=1)
        out.append(r)
    R = np.concatenate(out, axis=1)
    return np.where(np.isfinite(R), R, 1e3)
