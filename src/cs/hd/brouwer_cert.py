# HD refutation for every k: the "doubling chain" family. a_1 = 1, a_j ~ 2^(j-1) (2 <= j <= k-1), a_k = M (large).
# Members 2..k-1 Lucas-good: prod_{i != j}(1 - 2 a_i/a_j) = 1.  Fixed-point map T_j(a) = 2 a_{j-1} / (1 - 1/R_j(a)),
# R_j = prod_{i != j, i != j-1} (1 - 2 a_i / a_j).  Here: iterate T (numerics), then Krawczyk-certify the square system
# G_j(a) = prod_{i != j}(1 - 2 a_i/a_j) - 1 = 0 (j = 2..k-1) in unknowns a_2..a_{k-1}.  Members 1 and k checked dead.
import sys; sys.path.insert(0, '.')
from krk import certify
import mpmath
from mpmath import iv, mp
k = int(sys.argv[1]); M = int(sys.argv[2])
mp.dps = 60
def full(x): return [1] + list(x) + [M]
def G(x):
    a = full(x); out = []
    for j in range(1, k-1):
        p = 1
        for i in range(k):
            if i != j: p = p * (1 - 2 * a[i] / a[j])
        out.append(p - 1)
    return out
def J(x):
    a = full(x); n = k - 2; out = []
    for j in range(1, k-1):
        row = []
        for m in range(1, k-1):
            if m == j:
                # d/da_j prod_{i != j}(1 - 2 a_i/a_j) = sum_l (2 a_l / a_j^2) prod_{i != j, l}
                s = 0
                for l in range(k):
                    if l == j: continue
                    p = 2 * a[l] / a[j] ** 2
                    for i in range(k):
                        if i != j and i != l: p = p * (1 - 2 * a[i] / a[j])
                    s = s + p
                row.append(s)
            else:
                p = -2 / a[j]
                for i in range(k):
                    if i != j and i != m: p = p * (1 - 2 * a[i] / a[j])
                row.append(p)
        out.append(row)
    return out
# fixed-point iteration of T
x = [mpmath.mpf(2) ** (j - 1) for j in range(2, k)]
for it in range(200):
    a = full(x); new = []
    for j in range(1, k-1):
        R = mpmath.mpf(1)
        for i in range(k):
            if i != j and i != j - 1: R *= (1 - 2 * a[i] / a[j])
        new.append(2 * a[j-1] / (1 - 1 / R))
    x = new
mp.dps = int(sys.argv[3]) if len(sys.argv) > 3 else 60
ok, xs, K, X = certify(G, J, x, mpmath.mpf("1e-35"), dps=mp.dps)
a = [iv.mpf(1)] + X + [iv.mpf(M)]
order = all(a[i].b < a[i+1].a for i in range(k-1))
dead = []
for j in range(k):
    p = iv.mpf(1)
    for i in range(k):
        if i != j: p = p * (1 - 2 * a[i] / a[j])
    dead.append(not (p.a <= 1 <= p.b))
print(f"k={k} M={M}: Krawczyk={ok} ordered/distinct/positive={order} dead members (certified !=1): {[j+1 for j in range(k) if dead[j]]}")
print("   a =", [1] + [mpmath.nstr(v, 25) for v in xs] + [M])
print("   a_j / 2^(j-1) - 1 =", [mpmath.nstr(xs[j-1] / 2 ** j - 1, 6) for j in range(1, k-1)])
