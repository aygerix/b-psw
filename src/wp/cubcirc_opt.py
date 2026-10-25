"""Adversarial test of the CONJECTURE (real cubic normal form f = 2t+1+t(t+1)(a t + c), special points non-real):
for the member j with the smallest a (largest radius r = sqrt(2/a)), at least one of its upper special points
w_j ("-") and v_j ("+") is uncovered, i.e. max( min_{i != j} |f_i(w_j)|, min_{i != j} |f_i(v_j)| ) >= 1.
A counterexample needs only 3 members (one coverer per point). We MINIMISE
   M = max(|f_i(w_j)|, |f_i'(v_j)|)  over (a_j, c_j, a_i, c_i, a_i', c_i') with a_i, a_i' >= a_j,
all six special points non-real, by random restarts + Nelder-Mead (stdlib). Reports the smallest M found.
Usage: python3 cubcirc_opt.py restarts [seed]
"""
import sys, random, math

def f(a, c, t): return 2 * t + 1 + t * (t + 1) * (a * t + c)

def specials(a, c):
    d1 = (c - a) ** 2 - 8 * a; d2 = c * c - 8 * a
    if a <= 0 or d1 >= 0 or d2 >= 0: return None
    w = complex(-(c - a) / (2 * a), math.sqrt(-d1) / (2 * a)) - 1
    v = complex(-c / (2 * a), math.sqrt(-d2) / (2 * a))
    return w, v

def obj(x):
    aj, cj = math.exp(x[0]), x[1]
    ai, ci = aj * (1 + math.exp(x[2])), x[3]          # a_i > a_j
    ak, ck = aj * (1 + math.exp(x[4])), x[5]          # a_i' > a_j
    sj, si, sk = specials(aj, cj), specials(ai, ci), specials(ak, ck)
    if not (sj and si and sk): return 10.0
    w, v = sj
    return max(abs(f(ai, ci, w)), abs(f(ak, ck, v)))

def nelder_mead(fun, x0, step=0.5, iters=3000, tol=1e-13):
    n = len(x0); pts = [x0[:]] + [[x0[j] + (step if j == i else 0) for j in range(n)] for i in range(n)]
    vals = [fun(p) for p in pts]
    for _ in range(iters):
        order = sorted(range(n + 1), key=lambda i: vals[i]); pts = [pts[i] for i in order]; vals = [vals[i] for i in order]
        if vals[-1] - vals[0] < tol: break
        cen = [sum(p[j] for p in pts[:-1]) / n for j in range(n)]
        xr = [cen[j] + (cen[j] - pts[-1][j]) for j in range(n)]; fr = fun(xr)
        if fr < vals[0]:
            xe = [cen[j] + 2 * (cen[j] - pts[-1][j]) for j in range(n)]; fe = fun(xe)
            pts[-1], vals[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < vals[-2]:
            pts[-1], vals[-1] = xr, fr
        else:
            xc = [cen[j] + 0.5 * (pts[-1][j] - cen[j]) for j in range(n)]; fc = fun(xc)
            if fc < vals[-1]: pts[-1], vals[-1] = xc, fc
            else:
                for i in range(1, n + 1):
                    pts[i] = [pts[0][j] + 0.5 * (pts[i][j] - pts[0][j]) for j in range(n)]; vals[i] = fun(pts[i])
    i = min(range(n + 1), key=lambda i: vals[i]); return pts[i], vals[i]

def main():
    R = int(sys.argv[1]); rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    best = []
    for r in range(R):
        if r % 3 == 0:     # near-diagonal start: coverers close to member j
            aj = rng.uniform(-3, 3); cj = rng.uniform(-1, 1) * 2.5 * math.sqrt(8 * math.exp(aj)) + math.exp(aj) / 2
            x0 = [aj, cj, rng.uniform(-8, -3), cj + rng.gauss(0, 0.01), rng.uniform(-8, -3), cj + rng.gauss(0, 0.01)]
        else:
            x0 = [rng.uniform(-3, 3), rng.uniform(-6, 6), rng.uniform(-4, 3), rng.uniform(-6, 6), rng.uniform(-4, 3), rng.uniform(-6, 6)]
        if obj(x0) >= 10: continue
        x, v = nelder_mead(obj, x0)
        x, v = nelder_mead(obj, x, step=0.05)
        best.append((v, x))
    best.sort()
    print(f"{len(best)} feasible restarts; smallest M = max(|f_i(w_j)|, |f_i'(v_j)|) found:")
    for v, x in best[:8]:
        aj = math.exp(x[0]); print(f"  M = {v:.12f}   a_j = {aj:.5g}, c_j = {x[1]:.5g}, a_i/a_j = {1+math.exp(x[2]):.6g}, c_i = {x[3]:.6g}, a_i'/a_j = {1+math.exp(x[4]):.6g}, c_i' = {x[5]:.6g}")
    print("counterexamples (M < 1 - 1e-9):", sum(1 for v, _ in best if v < 1 - 1e-9))

if __name__ == "__main__":
    main()
