# Krawczyk certification of a real solution of a square polynomial system (mpmath interval arithmetic, outward rounding).
# Usage as a module: certify(G, J, x0, r) with G(x) -> list (works on mpf and mpi), J(x) -> matrix (list of lists).
import mpmath
from mpmath import iv, mp
def certify(G, J, x0, r, dps=60):
    mp.dps = dps; iv.dps = dps
    n = len(x0)
    x0 = [mpmath.mpf(v) for v in x0]
    # Newton polish in mp
    for _ in range(60):
        g = mpmath.matrix(G(x0)); Jm = mpmath.matrix(J(x0))
        dx = mpmath.lu_solve(Jm, g)
        x0 = [x0[i] - dx[i] for i in range(n)]
    Y = mpmath.inverse(mpmath.matrix(J(x0)))   # approximate inverse (point matrix)
    X = [iv.mpf([x0[i] - r, x0[i] + r]) for i in range(n)]
    xc = [iv.mpf(v) for v in x0]
    gc = G(xc)                       # interval enclosure of G(x0)
    JX = J(X)                        # interval Jacobian on box
    K = []
    for i in range(n):
        s = xc[i]
        for j in range(n): s = s - iv.mpf(Y[i, j]) * gc[j]
        for j in range(n):
            # (I - Y JX)_{ij}
            m = iv.mpf(1 if i == j else 0)
            for l in range(n): m = m - iv.mpf(Y[i, l]) * JX[l][j]
            s = s + m * (X[j] - xc[j])
        K.append(s)
    inside = all((K[i].a > X[i].a) and (K[i].b < X[i].b) for i in range(n))
    return inside, x0, K, X
