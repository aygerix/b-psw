"""Does the topmost special point of a real family escape every other member's sublevel set? (Theorem 10 extension test)
Members: real polynomials of even degree d with positive leading coefficient and f > 1 on R (no real special points),
f = c * prod_i ((t - x_i)^2 + b_i^2) with c * min_R prod > 1 enforced by scaling.
For random families of k members compute all special points (roots of f -+ 1), take the topmost one a* (member j*),
and test whether |f_i(a*)| < 1 for some i != j* ("covered"). Theorem 10 shows this never happens for d = 2.
Also report, for single members, whether the top of the upper sublevel set is a special point.
Usage: python3 topcover.py d k trials [seed]
"""
import sys, random, cmath, math

def pmul(p, q):
    r = [0j] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q): r[i + j] += a * b
    return r

def pval(p, t):
    v = 0j
    for c in reversed(p): v = v * t + c
    return v

def roots(p):                       # Durand-Kerner, p low->high, complex
    n = len(p) - 1; lc = p[-1]; q = [c / lc for c in p]
    z = [cmath.exp(2j * math.pi * (i + 0.25) / n) * (1 + abs(q[0]) ** (1 / n)) for i in range(n)]
    for _ in range(500):
        mx = 0
        for i in range(n):
            num = pval(q, z[i]); den = 1
            for j in range(n):
                if j != i: den *= (z[i] - z[j])
            dz = num / den; z[i] -= dz; mx = max(mx, abs(dz))
        if mx < 1e-14: break
    return z

def member(d, rng):
    p = [1 + 0j]
    for _ in range(d // 2):
        x = rng.uniform(-2, 2); b = math.exp(rng.uniform(-1.5, 1))
        p = pmul(p, [x * x + b * b, -2 * x, 1])
    mn = min(abs(pval(p, x / 200)) for x in range(-1000, 1001))
    s = math.exp(rng.uniform(0, 2.5)) / mn          # f > 1 on R (checked on a grid; min attained inside [-5,5])
    return [c * s for c in p]

def main():
    d, k, trials = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    rng = random.Random(int(sys.argv[4]) if len(sys.argv) > 4 else 1)
    covered = 0; examples = 0; nonspecial_top = 0; singles = 0
    for _ in range(trials):
        fam = [member(d, rng) for _ in range(k)]
        sp = []
        for j, f in enumerate(fam):
            for s in (-1, 1):
                g = f[:]; g[0] += s           # f + s: s = -1 -> f - 1 ("-" points), s = +1 -> f + 1 ("+" points)
                for r in roots(g):
                    if r.imag > 1e-9: sp.append((r.imag, j, r, s))
        h, js, a, s = max(sp)
        cov = [i for i, f in enumerate(fam) if i != js and abs(pval(f, a)) < 1]
        if cov:
            covered += 1
            if examples < 3:
                examples += 1
                print(f"COVERED: top special point {a:.4f} of member {js} ({'+' if s > 0 else '-'} point) inside sublevel set of members {cov}", flush=True)
        # single-member check: highest point of {|f|<=1} in the upper half plane vs highest special point of f
        f = fam[0]; best = 0.0
        for xi in range(-300, 301):
            x = xi / 100
            lo, hi = 0.0, 0.0
            for yi in range(1, 801):
                y = yi / 100
                if abs(pval(f, complex(x, y))) <= 1: hi = y
            best = max(best, hi)
        topsp = max(r.imag for s2 in (-1, 1) for r in roots([f[0] + s2] + f[1:]))
        singles += 1
        if best > topsp + 0.02: nonspecial_top += 1
    print(f"d = {d}, k = {k}: {trials} random families; topmost special point covered by another member in {covered}")
    print(f"single members whose sublevel set rises above their highest special point (grid 0.01): {nonspecial_top}/{singles}")

if __name__ == "__main__":
    main()
