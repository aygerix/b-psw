"""Real quadratic Williams families with a prescribed axis-class structure (Lemma 23). Real LM, stdlib only.

Members f(t) = alpha (t - x_c)^2 + beta with alpha > 0, beta > 1, grouped by axis x_c. Every special point is
x_c + i y with y > 0 (Lemma 23, step 1). A member is fixed by its two special heights y- < y+:
    alpha = 2 / (y+^2 - y-^2),   beta = 1 + alpha y-^2.
Unknowns: the axes x_2..x_(m-1) (x_1 = 0, x_m = 1) and log of each distinct height label in each class.
Residuals: F(x_c + i y) - 1 at "-" labels and F(x_c + i y) + 1 at "+" labels (complex, two real equations each),
plus a penalty when some member has y- >= y+, and box penalties (log-heights in [-6, 3], |axes| <= 3) that block the
degenerate limit heights -> infinity (all axes merge, alpha -> 0).
A zero residual would be a real quadratic Williams family with this structure.
Usage: python3 rq_struct.py STRUCT starts [seed]
  STRUCT: classes separated by '|', members separated by ',', each member "m<a>p<b>" (labels a, b local to the class);
  for example "m1p1|m1p1,m1p2|m1p1,m2p1|m1p1|m1p1" is [1 | 2- | 2+ | 1 | 1].
"""
import sys, random, math, cmath

def parse(s):
    classes = []
    for cs in s.split('|'):
        mem = []
        for ms in cs.split(','):
            a, b = ms[1:].split('p')
            mem.append((int(a), int(b)))
        classes.append(mem)
    return classes

def setup(classes):
    labels = []                                  # (class, sign, label id)
    for c, mem in enumerate(classes):
        for sgn in ('m', 'p'):
            ids = sorted({(a if sgn == 'm' else b) for a, b in mem})
            for i in ids: labels.append((c, sgn, i))
    return labels

def unpack(classes, labels, z):
    m = len(classes)
    xs = [0.0] + list(z[:m - 2]) + [1.0]
    h = {lab: math.exp(max(-30.0, min(30.0, v))) for lab, v in zip(labels, z[m - 2:])}
    mems, pen = [], 0.0
    for c, mem in enumerate(classes):
        for a, b in mem:
            ym, yp = h[(c, 'm', a)], h[(c, 'p', b)]
            d = yp * yp - ym * ym
            if d <= 1e-9:
                pen += 10.0 * (1e-9 - d + 1.0); d = 1e-9
            al = 2.0 / d
            mems.append((xs[c], al, 1.0 + al * ym * ym))
    return xs, h, mems, pen

def resid(classes, labels, z):
    xs, h, mems, pen = unpack(classes, labels, z)
    E = []
    for (c, sgn, i) in labels:
        t = complex(xs[c], h[(c, sgn, i)])
        F = 1 + 0j
        for x, al, be in mems: F *= al * (t - x) ** 2 + be
        target = 1.0 if sgn == 'm' else -1.0
        r = F - target
        E += [r.real, r.imag]
    if pen: E.append(pen)
    # keep away from the degenerate limit (all heights >> axis spread, alpha -> 0): confine log-heights to [-6, 3]
    m = len(classes)
    E.append(sum(max(0.0, v - 3.0) + max(0.0, -6.0 - v) for v in z[m - 2:]))
    E.append(sum(max(0.0, abs(v) - 3.0) for v in z[:m - 2]))
    return E

def solve(A, b):
    n = len(b); M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        p = max(range(c, n), key=lambda i: abs(M[i][c]))
        if abs(M[p][c]) < 1e-300: return None
        M[c], M[p] = M[p], M[c]
        for i in range(c + 1, n):
            f = M[i][c] / M[c][c]
            if f:
                for j in range(c, n + 1): M[i][j] -= f * M[c][j]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))) / M[i][i]
    return x

def lm(classes, labels, z, iters=300):
    n = len(z); lam = 1e-3
    E = resid(classes, labels, z); f = sum(e * e for e in E)
    for _ in range(iters):
        J = []
        for c in range(n):
            z2 = z[:]; hh = 1e-7 * max(1.0, abs(z[c])); z2[c] += hh
            E2 = resid(classes, labels, z2)
            m = min(len(E2), len(E))
            J.append([(E2[i] - E[i]) / hh for i in range(m)])
        m = min(len(r) for r in J)
        A = [[sum(J[a][i] * J[b][i] for i in range(m)) for b in range(n)] for a in range(n)]
        g = [-sum(J[a][i] * E[i] for i in range(m)) for a in range(n)]
        while True:
            A2 = [row[:] for row in A]
            for a in range(n): A2[a][a] += lam * (1 + A[a][a])
            dz = solve(A2, g)
            if dz is None: return z, f
            z2 = [a + b for a, b in zip(z, dz)]
            E2 = resid(classes, labels, z2); f2 = sum(e * e for e in E2)
            if f2 < f:
                z, E, f = z2, E2, f2; lam = max(lam / 3, 1e-15); break
            lam *= 4
            if lam > 1e10: return z, f
        if f < 1e-26: break
    return z, f

def main():
    classes = parse(sys.argv[1]); starts = int(sys.argv[2])
    rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    labels = setup(classes); m = len(classes)
    nunk = m - 2 + len(labels); neq = 2 * len(labels)
    print(f"structure {sys.argv[1]}: {m} classes, {len(labels)} heights, {nunk} unknowns, {neq} real equations", flush=True)
    best = []
    for s in range(starts):
        xs = sorted(rng.uniform(0.05, 0.95) for _ in range(m - 2))
        z = xs + [rng.uniform(-2.5, 1.5) for _ in labels]
        z, f = lm(classes, labels, z)
        best.append((f, z))
    best.sort(key=lambda u: u[0])
    print(f"smallest residual^2 over {starts} starts: " + ", ".join(f"{u[0]:.3e}" for u in best[:6]))
    f, z = best[0]
    xs, h, mems, pen = unpack(classes, labels, z)
    print("  best axes:", [round(v, 5) for v in xs])
    print("  best members (axis, alpha, beta):", [(round(a, 4), round(b, 5), round(c, 5)) for a, b, c in mems])
    print(f"  zeros (residual^2 < 1e-20): {sum(1 for u in best if u[0] < 1e-20)}")

if __name__ == "__main__":
    main()
