"""V2: does a free F/L side assignment beat Pomerance's mod-4 rule by more than a constant factor?

Input: freegen.c output (p, with A = (p-1)/2 and B = (p+1)/4 factored over the odd primes q <= T, q != 5, and no
side rule). For each p we compute o = ord_p 2 and z = z(p), keep class (1,2) (v2(o) = 1, v2(z) = 2), and record
the requirements: F(p) = prime powers of odd(o) (need n == 1 mod them), L(p) = prime powers of odd(z) (n == -1).
A side assignment s: {q} -> {F, L} admits p iff every q in F(p) has s(q) = F and every q in L(p) has s(q) = L.
For the admitted rows the subpool is optimised exactly by the min cut of src/pom/mincut.py (unit condition first).
Objective: the max ratio N / b(S) (Dinkelbach), and the lambda-optimum N - lambda b for a fixed lambda.
Search: start from the mod-4 rule (q == 1 mod 4 on F) and from random assignments; flip single primes while the
objective improves (first improvement, random order).
Usage: python3 sides.py freegen_output [restarts] [seed] [lambda]
"""
import sys, os, math, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pom"))
from mincut import solve, exact_b, clean_units
from profile2 import order_mod, zrank, v2

def parse(fn):
    rows = []
    for line in open(fn):
        if "A:" not in line: continue
        ps, rest = line.split(" A:")
        a_part, b_part = rest.split("B:")
        p = int(ps)
        if p < 7: continue
        fa = {int(t.split("^")[0]): int(t.split("^")[1]) for t in a_part.split()}
        fb = {int(t.split("^")[0]): int(t.split("^")[1]) for t in b_part.split()}
        fm = dict(fa); fm[2] = 1                      # p - 1 = 2 A
        fp = dict(fb); fp[2] = 2                      # p + 1 = 4 B
        o, z = order_mod(2, p, fm), zrank(p, fp)
        if (v2(o), v2(z)) != (1, 2): continue
        oo, zz = o >> 1, z >> 2
        F = {q: e for q, e in ((q, _val(oo, q)) for q in fa) if e}
        L = {q: e for q, e in ((q, _val(zz, q)) for q in fb) if e}
        rows.append((p, F, L))
    return rows

def _val(x, q):
    e = 0
    while x % q == 0: x //= q; e += 1
    return e

def admitted(rows, side):
    out = []
    for p, F, L in rows:
        if all(side[q] == 0 for q in F) and all(side[q] == 1 for q in L):
            out.append((p, {**F, **L}))
    return clean_units(out)

def best_ratio(sub):
    if not sub: return 0.0, 0, 0.0
    b0, _ = exact_b(sub)
    lam = len(sub) / max(b0, 1e-9); best = (0.0, 0, 0.0)
    for _ in range(30):
        S = solve(sub, lam); b, _ = exact_b(S)
        if not S or b < 1e-9: break
        r = len(S) / b
        if r > best[0] + 1e-12: best = (r, len(S), b)
        if abs(r - lam) < 1e-9: break
        lam = r
    return best

def lam_opt(sub, lam):
    if not sub: return 0.0, 0, 0.0
    S = solve(sub, lam); b, _ = exact_b(S)
    return len(S) - lam * b, len(S), b

def main():
    rows = parse(sys.argv[1])
    restarts = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    lam = float(sys.argv[4]) if len(sys.argv) > 4 else None
    qs = sorted({q for _, F, L in rows for q in list(F) + list(L)})
    print(f"class (1,2) rows: {len(rows)}; primes with requirements: {len(qs)}")
    nreq = [len(F) + len(L) for _, F, L in rows]
    print(f"requirements per row: mean {sum(nreq) / len(nreq):.2f}, max {max(nreq)}; "
          f"rows with <= 4: {sum(1 for x in nreq if x <= 4)}")
    def score(side):
        sub = admitted(rows, side)
        return (lam_opt(sub, lam) if lam is not None else best_ratio(sub)), len(sub)
    mod4 = {q: (0 if q % 4 == 1 else 1) for q in qs}
    (v, N, b), na = score(mod4)
    print(f"mod-4 rule: admitted {na}; objective {v:.4f} at N = {N}, b = {b:.2f}", flush=True)
    starts = [("mod4", dict(mod4))] + [(f"rand{r}", {q: rng.randrange(2) for q in qs}) for r in range(restarts)]
    overall = None
    for name, side in starts:
        (cur, N, b), na = score(side)
        improved = True
        while improved:
            improved = False
            order = qs[:]; rng.shuffle(order)
            for q in order:
                side[q] ^= 1
                (v, N2, b2), na2 = score(side)
                if v > cur + 1e-9:
                    cur, N, b, na = v, N2, b2, na2; improved = True
                else:
                    side[q] ^= 1
        nF = sum(1 for q in qs if side[q] == 0)
        agree = sum(1 for q in qs if side[q] == mod4[q])
        print(f"start {name}: local optimum objective {cur:.4f} at N = {N}, b = {b:.2f}; admitted {na}; "
              f"|F| = {nF}, agrees with mod-4 on {agree}/{len(qs)}", flush=True)
        if overall is None or cur > overall[0]: overall = (cur, N, b, dict(side))
    cur, N, b, side = overall
    print(f"BEST: objective {cur:.4f} at N = {N}, b = {b:.2f}; F side = {sorted(q for q in qs if side[q] == 0)}")

if __name__ == "__main__":
    main()
