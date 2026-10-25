# HD refutation, certificate 1: k = 5 linear common-root Korselt family f_i = 1 + a_i x with ALL a_i > 0 and exactly
# the three middle members Lucas-good. Fix a1 = 1 and a5 = 135 (rational); unknowns a2, a3, a4.
# Good_j: E_j = prod_{i != j}(a_j - 2 a_i) - a_j^4 = 0  (j = 2, 3, 4).
import sys; sys.path.insert(0, '.')
from krk import certify
import mpmath
from mpmath import iv
A1 = 1; A5 = int(sys.argv[1]) if len(sys.argv) > 1 else 135
def full(x): return [A1, x[0], x[1], x[2], A5]
def E(a, j):
    p = 1
    for i in range(5):
        if i != j: p = p * (a[j] - 2 * a[i])
    return p - a[j] ** 4
def G(x):
    a = full(x); return [E(a, 1), E(a, 2), E(a, 3)]
def J(x):
    # analytic partial derivatives of E_j w.r.t. a_m (m = 1, 2, 3 indices in a)
    a = full(x); out = []
    for j in (1, 2, 3):
        row = []
        for m in (1, 2, 3):
            if m == j:
                # d/da_j prod_{i != j}(a_j - 2 a_i) - 4 a_j^3
                s = 0
                for l in range(5):
                    if l == j: continue
                    p = 1
                    for i in range(5):
                        if i != j and i != l: p = p * (a[j] - 2 * a[i])
                    s = s + p
                row.append(s - 4 * a[j] ** 3)
            else:
                p = 1
                for i in range(5):
                    if i != j and i != m: p = p * (a[j] - 2 * a[i])
                row.append(-2 * p)
        out.append(row)
    return out
x0 = [float(v) for v in sys.argv[2:5]] if len(sys.argv) > 4 else [1.99924097798, 4.043784805681, 7.462060465493]
ok, xs, K, X = certify(G, J, x0, mpmath.mpf('1e-30'))
print("A5 =", A5, " Krawczyk K(X) in int(X):", ok)
print("a =", [A1] + [mpmath.nstr(v, 40) for v in xs] + [A5])
print("box radius 1e-30; box:", [str(v) for v in X])
# certified facts on the box: positivity, ordering, distinctness, and the dead members' residuals nonzero
a = [iv.mpf(A1)] + X + [iv.mpf(A5)]
print("ordering a1 < a2 < a3 < a4 < a5 on box:", all(a[i].b < a[i+1].a for i in range(4)))
for j in range(5):
    p = iv.mpf(1)
    for i in range(5):
        if i != j: p = p * (1 - 2 * a[i] / a[j])
    print(f"  member {j+1}: prod_(i!=j)(1 - 2 a_i/a_j) in {p}  -> {'GOOD (=1)' if p.a <= 1 <= p.b else 'DEAD (!=1)'}")
