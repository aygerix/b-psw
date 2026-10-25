"""Real cubic normal form f = 2t + 1 + t(t+1)(a t + c): the "-" points lie on |t + 1| = r, the "+" points on |t| = r,
r = sqrt(2/a) (root products of a s^2 + (c - a) s + 2 and a v^2 + c v + 2). Extremal test: in random families, is a special
point of the LARGEST-radius member ever covered (|f_i| < 1 for some other member i)? Also report for comparison the
topmost-point rule and the smallest-radius rule. Usage: python3 cubcirc.py trials [seed]
"""
import sys, random, math

def f(a, c, t): return 2 * t + 1 + t * (t + 1) * (a * t + c)

def specials(a, c):
    d1 = (c - a) ** 2 - 8 * a; d2 = c * c - 8 * a
    if d1 >= 0 or d2 >= 0: return None
    w = complex(-(c - a) / (2 * a), math.sqrt(-d1) / (2 * a)) - 1
    v = complex(-c / (2 * a), math.sqrt(-d2) / (2 * a))
    return w, v

def sample(rng):
    while True:
        a = math.exp(rng.uniform(-4, 3)); c = rng.uniform(-1, 1) * 3 * math.sqrt(8 * a) + a / 2
        sp = specials(a, c)
        if sp: return a, c, sp

def covered(fam, j, z):
    return any(abs(f(a, c, z)) < 1 for i, (a, c, sp) in enumerate(fam) if i != j)

def main():
    trials = int(sys.argv[1]); rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    cnt = {'largest-r': 0, 'smallest-r': 0, 'topmost': 0}; both = {'largest-r': 0, 'smallest-r': 0}
    for _ in range(trials):
        k = rng.choice([3, 5, 7, 9])
        fam = [sample(rng) for _ in range(k)]
        jl = min(range(k), key=lambda j: fam[j][0]); js = max(range(k), key=lambda j: fam[j][0])
        for name, j in (('largest-r', jl), ('smallest-r', js)):
            w, v = fam[j][2]
            cw, cv = covered(fam, j, w), covered(fam, j, v)
            if cw or cv: cnt[name] += 1
            if cw and cv: both[name] += 1
        pts = [(z.imag, j, z) for j, (a, c, sp) in enumerate(fam) for z in sp]
        h, jt, z = max(pts)
        if covered(fam, jt, z): cnt['topmost'] += 1
    print(f"{trials} random families (k in 3..9): families where the chosen member has a covered special point:")
    for kname in cnt: print(f"  {kname}: {cnt[kname]}" + (f"  (both its - and + points covered: {both[kname]})" if kname in both else ""))

if __name__ == "__main__":
    main()
