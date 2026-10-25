"""Symbolic D-basis data loader (pc7_cert_dsym.txt) -> fast numeric evaluation at float shapes, and full 7-member evaluation."""
import numpy as np, sympy as sp, itertools
e1, e2, lam, mu = sp.symbols('e1 e2 lam mu')
LAM3 = np.array([1, 0, 0, 0, -1, 0, 0, 0], float)
def _split(s):
    s = s.strip()[1:-1]; parts, dep, cur = [], 0, ''
    for ch in s:
        if ch in '([': dep += 1
        if ch in ')]': dep -= 1
        if ch == ',' and dep == 0: parts.append(cur); cur = ''
        else: cur += ch
    parts.append(cur); return parts
def load(path='pc7_cert_dsym.txt'):
    K = None; A = []; C = []
    for l in open(path):
        tag, rest = l.strip().split(' ', 1)
        ex = [sp.sympify(x.replace('^', '**'), locals=dict(e1=e1, e2=e2, lam=lam, mu=mu)) for x in _split(rest)]
        if tag == 'K': K = ex
        elif tag == 'A': A.append(ex)
        else: C.append(ex)
    allx = K + [x for r in A for x in r] + [x for r in C for x in r]
    f = sp.lambdify((e1, e2, lam, mu), allx, 'numpy')
    def at(r1, r2):
        E1, E2 = r1 + r2, r1 * r2; Dl = (1 - r1**2) * (1 - r2**2)
        v = np.array(f(E1, E2, 2 * E1 / Dl, 2 * (1 - E1**2 + E2) / Dl), float)
        return v[:8], v[8:80].reshape(8, 9), v[80:120].reshape(8, 5)
    return at
def Evec(K, A, C, ws):
    ws = np.asarray(ws, float); p = np.array([np.sum(ws**m) for m in range(0, 10)])
    pair = np.zeros(8)
    for n in range(2, 7):
        s = sum(0.5 * (p[a] * p[n - a] - p[n]) for a in range(1, n)); pair += C[:, n - 2] * s
    return K + A @ p[1:10] + pair + LAM3 * (p[1]**3 - 3 * p[1] * p[2] + 2 * p[3]) / 6
