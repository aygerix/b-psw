"""Pencil families f_a = L + a Q4 (real quartics with real special points rho1 < r1' < r2' < rho2, Lemma 25 draft).
L = cubic with L(rho1) = L(rho2) = 1, L(r1') = L(r2') = -1; Q4 = (t-rho1)(t-rho2)(t-r1')(t-r2').
For a > 0 admissible (g_a = (f_a - 1)/((t-rho1)(t-rho2)) and h_a = (f_a + 1)/((t-r1')(t-r2')) have non-real roots), record
sign Re Q4(omega) at the roots omega of g_a ("-" points) and Re Q4(nu) at the roots of h_a ("+" points).
Half-plane argument: the a_max member needs Re Q4(omega) > 0, Re Q4(nu) < 0; the a_min member needs the opposite.
If the signs are constant on the admissible a-range, no pencil family exists (any k).
Usage: python3 pencil_sign.py trials [seed]
"""
import sys, random, cmath

def polymul(p, q):
    r = [0.0] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q): r[i + j] += x * y
    return r

def polyval(p, x):          # p[0] + p[1] x + ...
    v = 0
    for c in reversed(p): v = v * x + c
    return v

def interp(xs, ys):          # Lagrange, coefficients low->high
    n = len(xs); res = [0.0] * n
    for i in range(n):
        num = [1.0]; den = 1.0
        for j in range(n):
            if j == i: continue
            num = polymul(num, [-xs[j], 1.0]); den *= xs[i] - xs[j]
        for d in range(n): res[d] += ys[i] * num[d] / den
    return res

def quad_roots(c0, c1, c2):
    disc = c1 * c1 - 4 * c2 * c0
    s = cmath.sqrt(disc)
    return (-c1 + s) / (2 * c2), (-c1 - s) / (2 * c2), disc

def divide_by(p, r1, r2):    # divide polynomial p (low->high) by (t-r1)(t-r2), return quotient (low->high)
    q = list(reversed(p))
    for r in (r1, r2):
        out = [q[0]]
        for c in q[1:]: out.append(c + out[-1] * r)
        q = out[:-1]
    return list(reversed(q))

def main():
    trials = int(sys.argv[1]); rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    mixed = 0; cfgs = 0
    for _ in range(trials):
        pts = sorted(rng.uniform(-3, 3) for _ in range(4))
        rho1, r1p, r2p, rho2 = pts
        L = interp([rho1, r1p, r2p, rho2], [1.0, -1.0, -1.0, 1.0])
        Q4 = polymul(polymul([-rho1, 1.0], [-rho2, 1.0]), polymul([-r1p, 1.0], [-r2p, 1.0]))
        signs_m, signs_p = set(), set(); adm = 0
        for i in range(1, 400):
            a = 10 ** (-3 + 6 * i / 400.0)
            f = [L[d] + a * Q4[d] if d < len(L) else a * Q4[d] for d in range(5)]
            fm = f[:]; fm[0] -= 1; fp = f[:]; fp[0] += 1
            g = divide_by(fm, rho1, rho2); h = divide_by(fp, r1p, r2p)
            w1, w2, dg = quad_roots(g[0], g[1], g[2]); v1, v2, dh = quad_roots(h[0], h[1], h[2])
            if dg >= 0 or dh >= 0: continue
            if g[2] <= 0 or h[2] <= 0: continue
            adm += 1
            signs_m.add(1 if polyval(Q4, w1).real > 0 else -1)
            signs_p.add(1 if polyval(Q4, v1).real > 0 else -1)
        if adm >= 2:
            cfgs += 1
            if len(signs_m) > 1 or len(signs_p) > 1:
                mixed += 1
                if mixed <= 5: print("MIXED signs:", [round(x, 3) for x in pts], "minus-signs", signs_m, "plus-signs", signs_p, "admissible a values", adm)
            elif cfgs <= 5:
                print("constant:", [round(x, 3) for x in pts], "Re Q4 at - points:", signs_m, "at + points:", signs_p, "admissible a values", adm)
    print(f"configurations with >= 2 admissible a: {cfgs}; with sign changes along a: {mixed}")

if __name__ == "__main__":
    main()
