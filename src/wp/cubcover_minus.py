"""Can EVERY member's upper "-" point be covered simultaneously (real cubic normal form)? Minimise
  M = max_j min_{i != j} |f_i(w_j)|  over k distinct valid members (a > 0, special points non-real).
M < 1 means the covering lemma alone cannot rule out the "-" half-system. Nelder-Mead from random starts (stdlib).
Usage: python3 cubcover_minus.py k restarts [seed]
"""
import sys, random, math
from cubcirc_opt import nelder_mead

def f(a, c, t): return 2 * t + 1 + t * (t + 1) * (a * t + c)

def wpoint(a, c):
    d1 = (c - a) ** 2 - 8 * a; d2 = c * c - 8 * a
    if a <= 0 or d1 >= 0 or d2 >= 0: return None
    return complex(-(c - a) / (2 * a), math.sqrt(-d1) / (2 * a)) - 1

def make(k):
    def obj(x):
        ms = [(math.exp(x[2 * j]), x[2 * j + 1]) for j in range(k)]
        ws = [wpoint(a, c) for a, c in ms]
        if any(w is None for w in ws): return 10.0
        sep = min(abs(ms[i][0] - ms[j][0]) + abs(ms[i][1] - ms[j][1]) for i in range(k) for j in range(i))
        if sep < 1e-6: return 10.0
        return max(min(abs(f(ms[i][0], ms[i][1], ws[j])) for i in range(k) if i != j) for j in range(k))
    return obj

def main():
    k, R = int(sys.argv[1]), int(sys.argv[2]); rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    obj = make(k); best = []
    for r in range(R):
        x0 = []
        for j in range(k): x0 += [rng.uniform(-6, 2), rng.uniform(-3, 3)]
        if obj(x0) >= 10: continue
        x, v = nelder_mead(obj, x0, step=0.5, iters=6000); x, v = nelder_mead(obj, x, step=0.05, iters=6000)
        best.append((v, x))
    best.sort()
    print(f"k = {k}: {len(best)} feasible restarts; smallest M = max_j min_i |f_i(w_j)|:")
    for v, x in best[:5]:
        print(f"  M = {v:.6f}  (a, c) = {[(round(math.exp(x[2*j]), 6), round(x[2*j+1], 6)) for j in range(k)]}")
    print("all '-' points covered (M < 1):", sum(1 for v, _ in best if v < 1 - 1e-9))

if __name__ == "__main__":
    main()
