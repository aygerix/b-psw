"""Shared loader for pc7_cert_pw.txt (from pc7_cert_pw.gp): pointwise boundary functionals E1..E8 of the k = 7 (2,2) pencil
in 3-body form E_i = kap_i + sum_j a_i(w_j) + sum_{i<j} b_i(w_i, w_j) + lam3_i e3(w)."""
import re, numpy as np, sympy as sp
r1, r2, w, v, e1, e2, lam, mu = sp.symbols('r1 r2 w v e1 e2 lam mu')
Dl = (1 - r1**2) * (1 - r2**2)
SUBS = {e1: r1 + r2, e2: r1 * r2, lam: 2 * (r1 + r2) / Dl, mu: 2 * (1 - (r1 + r2)**2 + r1 * r2) / Dl}
LAM3 = np.array([1, 1, 0, 0, -1, -1, 0, 0], float)
NAMES = ["dV(r1)", "dV(r2)", "dV'(r1)", "dV'(r2)", "dW(1)", "dW(-1)", "dW'(1)", "dW'(-1)"]
def _sym(s, sub=True):
    e = sp.sympify(s.replace('^', '**'), locals={'r1': r1, 'r2': r2, 'w': w, 'v': v, 'e1': e1, 'e2': e2, 'lam': lam, 'mu': mu})
    return e.subs(SUBS) if sub else e
def load(path='pc7_cert_pw.txt'):
    lines = [l.strip() for l in open(path) if l.strip()]
    d = {}
    for l in lines:
        k, s = l.split(' = ', 1); d[k] = s
    # robust top-level split
    def split_top(s):
        s = s.strip()[1:-1]; parts, dep, cur = [], 0, ''
        for ch in s:
            if ch in '([': dep += 1
            if ch in ')]': dep -= 1
            if ch == ',' and dep == 0: parts.append(cur); cur = ''
            else: cur += ch
        parts.append(cur); return parts
    kap = [_sym(x) for x in split_top(d['kap'])]
    A = []
    for i in range(1, 9):
        cs = [_sym(x) for x in split_top(d['a%d' % i])]       # highest first: w^n .. w^1
        n = len(cs); A.append([(n - m, cs[m]) for m in range(n)])  # (power, coeff)
    Bc = []
    for i in range(1, 9):
        P = sp.Poly(sp.expand(_sym(d['b%d' % i]) * 1), w, v) if False else None
        Pn = sp.Poly(sp.expand(_sym(d['b%d' % i], sub=False)), w, v)
        Bc.append([(m, n, c.subs(SUBS)) for (m, n), c in Pn.terms()])
    f = lambda e: sp.lambdify((r1, r2), e, 'math')
    return dict(kap=[f(e) for e in kap], A=[[(p, f(c)) for p, c in a] for a in A], B=[[(m, n, f(c)) for m, n, c in b] for b in Bc],
                sym=dict(kap=kap, A=A, B=Bc))
def at_shape(D, R1, R2):
    """numeric data at a shape: kap (8,), Acoef (8, 10) [power 0..9], Bcoef (8, 10, 10)"""
    kap = np.array([f(R1, R2) for f in D['kap']])
    Ac = np.zeros((8, 10)); Bc = np.zeros((8, 10, 10))
    for i in range(8):
        for p, f in D['A'][i]: Ac[i, p] = f(R1, R2)
        for m, n, f in D['B'][i]: Bc[i, m, n] = f(R1, R2)
    return kap, Ac, Bc
