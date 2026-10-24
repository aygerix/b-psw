"""Does the top of the sublevel set {|f| <= 1} (upper half plane) of a real cubic in the Lemma 24 normal form
f = a t^3 + b t^2 + (b - a + 2) t + 1 (a > 0; g = (f-1)/t and h = (f+1)/(t+1) without real roots)
coincide with its highest non-real special point?  If yes for all such cubics, Theorem 10's modulus argument extends.
Grid test; also tests random k-member families for a covered topmost non-real special point.
Usage: python3 cubtop.py trials [seed]
"""
import sys, random, cmath, math

def f(a, b, t): return ((a * t + b) * t + (b - a + 2)) * t + 1

def specials(a, b):
    out = []
    for (A, B, C, s) in ((a, b, b - a + 2, -1), (a, b - a, 2.0, 1)):
        disc = B * B - 4 * A * C
        if disc >= 0: return None
        out.append((complex(-B / (2 * A), math.sqrt(-disc) / (2 * A)), s))
    return out

def sample(rng):
    while True:
        a = math.exp(rng.uniform(-3, 3)); b = rng.uniform(-10, 10)
        sp = specials(a, b)
        if sp: return a, b, sp

def top_of_sublevel(a, b, ymax):
    best = 0.0
    for xi in range(-400, 401):
        x = xi / 100.0
        for yi in range(1, 401):
            y = ymax * yi / 400.0
            if abs(f(a, b, complex(x, y))) <= 1: best = max(best, y)
    return best

def main():
    trials = int(sys.argv[1]); rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    above = 0; ex = 0
    for _ in range(trials):
        a, b, sp = sample(rng)
        hs = max(z.imag for z, s in sp)
        top = top_of_sublevel(a, b, 3 * hs + 1)
        if top > hs * 1.01 + 0.01:
            above += 1
            if ex < 4: ex += 1; print(f"a={a:.4g} b={b:.4g}: sublevel top {top:.4f} > highest special {hs:.4f}", flush=True)
    print(f"single cubics with sublevel set above their highest non-real special point: {above}/{trials}")
    cov = 0
    for _ in range(2000):
        k = rng.choice([3, 5, 7])
        fam = [sample(rng) for _ in range(k)]
        pts = [(z.imag, j, z) for j, (a, b, sp) in enumerate(fam) for z, s in sp]
        h, js, z = max(pts)
        if any(abs(f(a, b, z)) < 1 for j, (a, b, sp) in enumerate(fam) if j != js): cov += 1
    print(f"random families: topmost non-real special point covered by another member in {cov}/2000")

if __name__ == "__main__":
    main()
