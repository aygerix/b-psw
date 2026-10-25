"""Numeric check of the (2,2) quartic pencil reduction used in gen_quartic22.py: for random u1, u2 in (-1, 1) and random
c_i > 0, at the roots w of K-_j (resp. K+_j) the product F(w) equals prod_i (1 + (1 - c_i/c_j) X(w)) with X = f_0 - 1
(resp. prod_i (-1 + (1 - c_i/c_j) Y(w)), Y = f_0 + 1). Prints the max discrepancy."""
import random, cmath
def polymul(a, b):
    r = [0j] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): r[i + j] += x * y
    return r
def ev(p, t): return sum(c * t ** i for i, c in enumerate(p))
def roots2(q):  # q = [c0, c1, c2]
    c0, c1, c2 = q; d = cmath.sqrt(c1 * c1 - 4 * c2 * c0); return [(-c1 + d) / (2 * c2), (-c1 - d) / (2 * c2)]
rng = random.Random(5); worst = 0
for trial in range(200):
    u1, u2 = sorted(rng.uniform(-0.95, 0.95) for _ in range(2))
    G = [-1, 0, 1]; H = [u1 * u2, -(u1 + u2), 1]
    # f0 = 1 + G (al t + be) with f0(u) = -1 at u1, u2: (u^2 - 1)(al u + be) = -2
    a11, a12, b1 = (u1 * u1 - 1) * u1, (u1 * u1 - 1), -2
    a21, a22, b2 = (u2 * u2 - 1) * u2, (u2 * u2 - 1), -2
    det = a11 * a22 - a12 * a21; al = (b1 * a22 - a12 * b2) / det; be = (a11 * b2 - b1 * a21) / det
    ell = [be, al]; f0 = [x + (1 if i == 0 else 0) for i, x in enumerate(polymul(G, ell) + [0])][:4]
    P = polymul(G, H)
    k = 5; c = [rng.uniform(0.2, 3) for _ in range(k)]
    fs = [[f0[i] if i < 4 else 0 for i in range(5)] for _ in range(k)]
    for i in range(k): fs[i] = [fs[i][m] + c[i] * P[m] for m in range(5)]
    j = rng.randrange(k)
    Km = [ell[0] + c[j] * H[0], ell[1] + c[j] * H[1], c[j] * H[2]]
    # r0 = (f0 + 1) / H (linear): solve via two points
    for w in roots2(Km):
        F = 1
        for f in fs: F *= ev(f, w)
        X = ev(f0, w) - 1
        pred = 1
        for i in range(k): pred *= 1 + (1 - c[i] / c[j]) * X
        worst = max(worst, abs(F - pred) / (1 + abs(F)))
    # "+" side: K+_j = r0 + c_j G with f0 + 1 = H r0
    # r0 coefficients: f0 + 1 = (t^2 - sg t + rh)(ga t + de): ga = al, de from constant term: rh de = f0(0) + 1
    ga = al; de = (f0[0] + 1) / H[0]
    Kp = [de + c[j] * G[0], ga, c[j]]
    for v in roots2(Kp):
        F = 1
        for f in fs: F *= ev(f, v)
        Y = ev(f0, v) + 1
        pred = 1
        for i in range(k): pred *= -1 + (1 - c[i] / c[j]) * Y
        worst = max(worst, abs(F - pred) / (1 + abs(F)))
print("max relative discrepancy over 200 random pencils:", worst)
