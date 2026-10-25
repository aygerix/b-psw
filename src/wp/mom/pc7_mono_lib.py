"""pc7_mono: structure of the k = 7 (2,2) pencil boundary expressions.
B_i = C0_i + sum_j phi1_i(w_j) + sum_{j<l} phi2_i(w_j, w_l) + kap_i e3(w), with
phi2_i(u, v) = u v sum_n A_{i,n} h_n(u, v)  (h_n complete homogeneous), i.e. phi2 = u v [Gam(u) - Gam(v)]/(u - v),
Gam_i(x) = sum_n A_{i,n} x^(n+1).  Table pc7_mono_tab.txt from pc7_mono_gen.gp."""
import numpy as np
from fractions import Fraction
TAB = []
for line in open('pc7_mono_tab.txt'):
    s = line.split(None, 4)
    if len(s) < 5: continue
    i, a, b, c = map(int, s[:4]); TAB.append((i, a, b, c, compile(s[4].strip().replace('^', '**'), 'c', 'eval')))
def shape_env(r1, r2, exact=False):
    if exact: r1, r2 = Fraction(r1), Fraction(r2)
    e1, e2 = r1 + r2, r1 * r2; Dl = (1 - r1**2) * (1 - r2**2)
    return dict(e1=e1, e2=e2, lam=2 * e1 / Dl, mu=2 * (1 - e1**2 + e2) / Dl)
def parts(r1, r2, exact=False):
    """returns C0[8], phi1[8][10] (coef of w^a, a = 0..9; a = 0 unused), A[8][5] (pair h_n coefs), kap[8]"""
    env = shape_env(r1, r2, exact)
    z = (lambda: Fraction(0)) if exact else (lambda: 0.0)
    C0 = [z() for _ in range(8)]; P1 = [[z() for _ in range(10)] for _ in range(8)]
    A = [[z() for _ in range(5)] for _ in range(8)]; kap = [z() for _ in range(8)]
    seen = {}
    for i, a, b, c, code in TAB:
        v = eval(code, {}, env); i -= 1
        if a == b == c == 0: C0[i] = v
        elif b == 0: P1[i][a] = v
        elif c == 0:
            n = a + b - 2   # u^a v^b with a >= b >= 1 belongs to u v h_(a+b-2)
            if (i, n) in seen: assert seen[(i, n)] == v, (i, a, b, n)
            seen[(i, n)] = v; A[i][n] = v
        else: kap[i] = v
    return C0, P1, A, kap
def evalB(r1, r2, w):
    C0, P1, A, kap = parts(r1, r2)
    w = np.asarray(w, float); out = []
    for i in range(8):
        s = C0[i] + sum(P1[i][a] * np.sum(w**a) for a in range(1, 10))
        for j in range(len(w)):
            for l in range(j + 1, len(w)):
                u, v = w[j], w[l]
                s += u * v * sum(A[i][n] * sum(u**m * v**(n - m) for m in range(n + 1)) for n in range(5))
        e3 = sum(w[a] * w[b] * w[c] for a in range(len(w)) for b in range(a + 1, len(w)) for c in range(b + 1, len(w)))
        out.append(s + kap[i] * e3)
    return np.array(out)
