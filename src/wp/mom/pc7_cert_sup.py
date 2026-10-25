"""Pair relaxation: which supports suffice, is lam3 = 0 possible, and what do L1-minimal alpha look like?
pi_alpha(w, v) = kap/21 + (a(w) + a(v))/6 + b(w, v) >= 0 on (0, inf)^2, kap_alpha = 1, lam3_alpha >= 0."""
import sys, time, itertools, numpy as np
from scipy.optimize import linprog
from pc7_cert_lib import load, at_shape, LAM3, NAMES
D = load()
PW = np.exp(np.linspace(np.log(1e-5), np.log(1e5), 46))
WV = np.array([(a, b) for a in PW for b in PW if b >= a])
def aval(Ac, x): return np.array([np.polyval(Ac[i][::-1], x) for i in range(8)])
def bval(Bc, a, b):
    M = np.outer(a ** np.arange(10), b ** np.arange(10)); return np.array([np.sum(Bc[i] * M) for i in range(8)])
def nrm(R): return R / np.maximum(np.abs(R).max(axis=1, keepdims=True), 1e-300)
cache = {}
def rowsP(R1, R2):
    if (R1, R2) not in cache:
        kap, Ac, Bc = at_shape(D, R1, R2)
        cache[(R1, R2)] = (kap, nrm(np.array([kap / 21 + (aval(Ac, a) + aval(Ac, b)) / 6 + bval(Bc, a, b) for a, b in WV])))
    return cache[(R1, R2)]
def lp(R1, R2, supp=range(8), lam0=False, l1=False):
    kap, R = rowsP(R1, R2); s = list(supp); n = len(s)
    if l1:   # alpha = p - q
        c = np.ones(2 * n); Aub = -np.hstack([R[:, s], -R[:, s]]); Aub = np.vstack([Aub, -np.hstack([LAM3[s], -LAM3[s]])])
        Aeq = [np.concatenate([kap[s], -kap[s]])]; beq = [1.0]
        if lam0: Aeq.append(np.concatenate([LAM3[s], -LAM3[s]])); beq.append(0.0)
        r = linprog(c, A_ub=Aub, b_ub=np.zeros(len(Aub)), A_eq=np.array(Aeq), b_eq=beq, bounds=[(0, 1e5)] * (2 * n), method='highs')
        if r.status != 0: return None
        a = np.zeros(8); a[s] = r.x[:n] - r.x[n:]; return a
    Aub = -np.vstack([R[:, s], LAM3[s][None, :]]); Aeq = [kap[s]]; beq = [1.0]
    if lam0: Aeq.append(LAM3[s]); beq.append(0.0)
    r = linprog(np.zeros(n), A_ub=Aub, b_ub=np.zeros(len(Aub)), A_eq=np.array(Aeq), b_eq=beq, bounds=[(-1e5, 1e5)] * n, method='highs')
    return r.status == 0
grid = np.linspace(-0.9, 0.9, 7)
shapes = [(a, b) for a in grid for b in grid if b >= a - 1e-12]
print("lam3 = 0 forced:", sum(lp(*s, lam0=True) for s in shapes), "/", len(shapes), flush=True)
for name, supp in [("no derivatives", [0, 1, 4, 5]), ("V only", [0, 1, 2, 3]), ("W only", [4, 5, 6, 7]), ("values + V'", [0, 1, 2, 3, 4, 5]), ("values + W'", [0, 1, 4, 5, 6, 7])]:
    print(f"support {name}:", sum(lp(*s, supp=supp) for s in shapes), "/", len(shapes), flush=True)
good = []
for size in (4, 5, 6):
    for supp in itertools.combinations(range(8), size):
        if all(lp(*s, supp=supp) for s in shapes): good.append(tuple(i + 1 for i in supp))
    print(f"supports of size {size} working at all shapes: {[g for g in good if len(g) == size]}", flush=True)
print("L1-minimal alpha (functionals:", NAMES, ")")
for s in [(-0.6, 0.6), (-0.6, 0.3), (-0.6, 0.0), (-0.3, -0.3), (0.0, 0.0), (0.3, 0.6), (-0.9, 0.9), (-0.9, 0.6)]:
    a = lp(*s, l1=True); print(s, None if a is None else np.round(a, 5))
