#!/usr/bin/env python3
"""xn line: extra necessary conditions for real cubic Williams families (core evaluation).

Normal form f_i = 2t + 1 + t(t+1)(a_i t + c_i), a_i > 0; optional linear member 2t + 1 ("lin").
g_i = (f_i - 1)/t, P = prod g_i, V := (F - 1)/(t P).  Moment variables v = 1/a, u = c/a.

Exact reformulation used here (elementary, see notes in xn_report):
  V is a polynomial  <=>  the Taylor coefficients at tau = t + 1 = 0 of
        E(tau) := Vpoly(tau) - (1 - tau)^(k-1) G(tau),   G = prod(1 + x_i) - prod x_i,  x_i = 1/(f_i - 1),
  vanish for orders 0..2k-1, where Vpoly = polynomial part of t^(k-1) prod(1 + x_i) at t = oo
  (E = -R/P with R = remainder of (F-1)/t mod P, deg R < 2k, P(-1) = 2^k != 0).
  Truncation "order N" = equations E_0..E_{N-1}; N = 2 is the four boundary equations (with the reflection).
  The reflected family (c -> a - c, i.e. u -> 1 - u) gives the h-side conditions h_i | (F + 1)/(t + 1).
All routines are generic in the number type (Fraction for exact checks, float for numerics).
"""
from fractions import Fraction


def smul(A, B, N):
    C = [A[0] * 0] * N
    for i, a in enumerate(A[:N]):
        if a == 0:
            continue
        for j, b in enumerate(B[:N - i]):
            C[i + j] += a * b
    return C


def sinv(A, N):
    """1/A mod tau^N, A[0] != 0"""
    one = A[0] / A[0]
    R = [A[0] * 0] * N
    R[0] = one / A[0]
    for n in range(1, N):
        acc = A[0] * 0
        for j in range(1, min(n, len(A) - 1) + 1):
            acc += A[j] * R[n - j]
        R[n] = -acc / A[0]
    return R


def xs(u, v, K):
    """coefficients 0..K-1 (in s = 1/t) of x = v s^3 / (1 + (1+u) s + (u + 2v) s^2)"""
    z = u * 0
    d = [z] * K
    if K > 0:
        d[0] = z + 1
    if K > 1:
        d[1] = -(1 + u)
    for j in range(2, K):
        d[j] = -(1 + u) * d[j - 1] - (u + 2 * v) * d[j - 2]
    c = [z] * K
    for n in range(3, K):
        c[n] = v * d[n - 3]
    return c


def vpoly_s(U, Vv, lin, k):
    """T_j, j = 0..k-1: Vpoly(t) = sum_j T_j t^(k-1-j) (T = prod(1 + x_i) mod s^k)"""
    z = U[0] * 0 if len(U) else Fraction(0)
    T = [z] * k
    T[0] = z + 1
    for u, v in zip(U, Vv):
        X = xs(u, v, k)
        X[0] += 1
        T = smul(T, X, k)
    if lin:
        L = [z] * k
        L[0] = z + 1
        if k > 1:
            L[1] = z + Fraction(1, 2) if isinstance(z, Fraction) else z + 0.5
        T = smul(T, L, k)
    return T


def binom(n, r):
    from math import comb
    return comb(n, r) if 0 <= r <= n else 0


def vpoly_tau(T, k):
    """coefficients in tau of sum_j T_j (tau - 1)^(k-1-j)"""
    z = T[0] * 0
    out = [z] * k
    for j, tj in enumerate(T):
        e = k - 1 - j
        for n in range(e + 1):
            out[n] += tj * binom(e, n) * (-1) ** (e - n)
    return out


def G_tau(A, C, lin, k, N):
    """Taylor coefficients 0..N-1 at tau = 0 of (1 - tau)^(k-1) * G(tau)"""
    z = A[0] * 0 if len(A) else Fraction(0)
    one = z + 1
    P1 = [z] * N
    P1[0] = one
    P2 = [z] * N
    P2[0] = one
    mem = [(a, c) for a, c in zip(A, C)]
    if lin:
        mem.append((z, z))
    for a, c in mem:
        fm1 = [z - 2, 2 + a - c, c - 2 * a, a] + [z] * max(0, N - 4)   # f - 1 in tau
        X = sinv(fm1[:max(N, 1)], N)
        X1 = list(X)
        X1[0] += 1
        P1 = smul(P1, X1, N)
        P2 = smul(P2, X, N)
    G = [P1[n] - P2[n] for n in range(N)]
    B = [one * binom(k - 1, n) * (-1) ** n for n in range(N)]
    return smul(B, G, N)


def E_eqs(U, Vv, lin, N):
    """E_0..E_{N-1} for the family with moment data (u_i, v_i)"""
    k = len(U) + lin
    A = [1 / v for v in Vv]
    C = [u / v for u, v in zip(U, Vv)]
    T = vpoly_s(U, Vv, lin, k)
    Vt = vpoly_tau(T, k)
    z = U[0] * 0
    Vt = Vt + [z] * max(0, N - k)
    Gt = G_tau(A, C, lin, k, N)
    return [Vt[n] - Gt[n] for n in range(N)]


def E_both(U, Vv, lin, N):
    return E_eqs(U, Vv, lin, N) + E_eqs([1 - u for u in U], Vv, lin, N)


# ---------------- exact cross-check (polynomial division) -----------------
def pmul(A, B):
    C = [A[0] * 0] * (len(A) + len(B) - 1)
    for i, a in enumerate(A):
        for j, b in enumerate(B):
            C[i + j] += a * b
    return C


def pdivmod(Nm, D):
    Nm = list(Nm)
    D = list(D)
    while len(D) > 1 and D[-1] == 0:
        D.pop()
    while len(Nm) > 1 and Nm[-1] == 0:
        Nm.pop()
    q = [Nm[0] * 0] * max(1, len(Nm) - len(D) + 1)
    while len(Nm) >= len(D) and any(x != 0 for x in Nm):
        d = len(Nm) - len(D)
        co = Nm[-1] / D[-1]
        q[d] = co
        for i, b in enumerate(D):
            Nm[i + d] -= co * b
        Nm.pop()
    return q, Nm


def exact_check(U, Vv, lin, N):
    """returns (E via series, -R/P series via division, quotient == Vpoly?)"""
    k = len(U) + lin
    A = [1 / v for v in Vv]
    C = [u / v for u, v in zip(U, Vv)]
    mem = [(a, c) for a, c in zip(A, C)] + ([(Fraction(0), Fraction(0))] if lin else [])
    F = [Fraction(1)]
    P = [Fraction(1)]
    for a, c in mem:
        f = [Fraction(-1), 2 + a - c, c - 2 * a, a]          # f in tau
        F = pmul(F, f)
        fm1 = [Fraction(-2), 2 + a - c, c - 2 * a, a]
        g, r = pdivmod(fm1, [Fraction(-1), Fraction(1)])      # divide by t = tau - 1
        assert all(x == 0 for x in r)
        P = pmul(P, g)
    F[0] -= 1
    Nn, r = pdivmod(F, [Fraction(-1), Fraction(1)])
    assert all(x == 0 for x in r)
    Q, R = pdivmod(Nn, P)
    R = R + [Fraction(0)] * (N + 1)
    ser = smul(R[:N], sinv(P[:N] + [Fraction(0)] * max(0, N - len(P)), N), N)
    T = vpoly_s(U, Vv, lin, k)
    Vt = vpoly_tau(T, k)
    return E_eqs(U, Vv, lin, N), [-x for x in ser], (Q + [Fraction(0)] * (k - len(Q)))[:k] == Vt[:k]


if __name__ == "__main__":
    import random
    random.seed(1)
    for (m, lin) in [(3, 0), (2, 1), (5, 0), (4, 1), (7, 0), (8, 1), (9, 0)]:
        U = [Fraction(random.randint(-20, 20), random.randint(1, 9)) for _ in range(m)]
        Vv = [Fraction(random.randint(1, 30), random.randint(1, 9)) for _ in range(m)]
        k = m + lin
        e1, e2, okq = exact_check(U, Vv, lin, 2 * k)
        print(m, lin, "series == -R/P:", e1 == e2, " quotient == Vpoly:", okq)
