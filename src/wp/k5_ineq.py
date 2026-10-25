"""Feasibility of the necessary conditions for a REAL cubic Williams family with k = 5 (normal form), derived from
the fixed Laurent cofactor V = y^4 + 15/4 y + 45/16 (NOTEBOOK session 8):
  s_i = psi_i'(-1) > 0, alpha_i = psi_i'(0) > 0, b_i = alpha_i + s_i - 4 > 0 (leading coefficients),
  (1) e2(s) = 360, (2) e1(s) = 45/16 e5(s), (3) e2(alpha) = 360, (4) e1(alpha) = 45/16 e5(alpha),
  (5) sum 1/b_i = 15/4,
  (6) non-real "-" points: (2 - s_i)^2 < 8 b_i;  (7) non-real "+" points: (2 - alpha_i)^2 < 8 b_i.
Minimises the total violation (log-parametrised, Nelder-Mead from random starts). A strictly positive minimum
over many starts is EVIDENCE of infeasibility (then try to prove it); a zero gives a point to feed the full system.
Usage: python3 k5_ineq.py starts [seed]
"""
import sys, math, random
from cubcirc_opt import nelder_mead

def esym(v):
    e = [1.0] + [0.0] * len(v)
    for x in v:
        for r in range(len(v), 0, -1): e[r] += e[r - 1] * x
    return e

def viol(z):
    s = [math.exp(t) for t in z[:5]]; al = [math.exp(t) for t in z[5:]]
    b = [al[i] + s[i] - 4 for i in range(5)]
    es, ea = esym(s), esym(al)
    v = 0.0
    v += (math.log(es[2] / 360)) ** 2 + (math.log(es[1] / (45 / 16 * es[5]))) ** 2
    v += (math.log(ea[2] / 360)) ** 2 + (math.log(ea[1] / (45 / 16 * ea[5]))) ** 2
    for i in range(5):
        if b[i] <= 1e-9: v += (1 + (1e-9 - b[i])) ** 2 * 10; continue
        v += max(0.0, (2 - s[i]) ** 2 / (8 * b[i]) - 1 + 1e-6) ** 2
        v += max(0.0, (2 - al[i]) ** 2 / (8 * b[i]) - 1 + 1e-6) ** 2
    if all(x > 0 for x in b):
        v += (math.log(sum(1 / x for x in b) / 3.75)) ** 2
    return v

def main():
    starts = int(sys.argv[1]); rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    best = []
    for st in range(starts):
        z0 = [rng.gauss(math.log(4), 1.5) for _ in range(10)]
        z, f = nelder_mead(viol, z0, step=0.7, iters=6000)
        z, f = nelder_mead(viol, z, step=0.1, iters=6000)
        best.append((f, z))
    best.sort(key=lambda t: t[0])
    print(f"{starts} starts; smallest violations:", [round(f, 8) for f, _ in best[:8]])
    f, z = best[0]
    s = [math.exp(t) for t in z[:5]]; al = [math.exp(t) for t in z[5:]]
    print("best s =", [round(x, 5) for x in s]); print("best alpha =", [round(x, 5) for x in al])
    print("b =", [round(al[i] + s[i] - 4, 5) for i in range(5)])

if __name__ == "__main__":
    main()
