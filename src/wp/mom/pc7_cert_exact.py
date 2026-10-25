"""Exact (Fraction) per-shape data from pc7_cert_shape.gp dumps; rows of the pair / full / uni relaxations computed exactly
at rational grid points, then normalised and passed to the LP as floats."""
import sys, itertools, numpy as np
from fractions import Fraction as Fr
from scipy.optimize import linprog
LAM3 = [1, 0, 0, 0, -1, 0, 0, 0]   # B- and D-basis (pc7_cert_shape.gp BASIS = 1, 2)
def parse_vec(s):
    s = s.strip().strip('[]'); return [Fr(x.strip()) for x in s.split(',')]
def load(path):
    out = []; cur = None
    for l in open(path):
        l = l.strip()
        if not l: continue
        tag, rest = l.split(' ', 1)
        if tag == 'S': cur = dict(r=tuple(Fr(x) for x in rest.split()), K=None, A=[], C=[]); out.append(cur)
        elif tag == 'K': cur['K'] = parse_vec(rest)
        elif tag == 'A': cur['A'].append(parse_vec(rest))
        elif tag == 'C': cur['C'].append(parse_vec(rest))
    return out
def aval(S, i, x): return sum(S['A'][i][m - 1] * x**m for m in range(1, 10))
def bval(S, i, x, y): return sum(S['C'][i][n - 2] * sum(x**a * y**(n - a) for a in range(1, n)) for n in range(2, 7))
def grid(lo=-12, hi=12, step=1):  # powers of 2 times 1, 3/2
    g = []
    for e in range(lo, hi + 1, step): g += [Fr(2)**e, Fr(3, 2) * Fr(2)**e]
    return sorted(set(g))
def nrm(rows):
    R = np.array([[float(x) for x in r] for r in rows]); return R / np.maximum(np.abs(R).max(axis=1, keepdims=True), 1e-300)
def pair_rows(S, G):
    rows = []
    for x, y in itertools.combinations_with_replacement(G, 2):
        ax = [aval(S, i, x) for i in range(8)]; ay = [aval(S, i, y) for i in range(8)]
        rows.append([S['K'][i] / 21 + (ax[i] + ay[i]) / 6 + bval(S, i, x, y) for i in range(8)])
    return rows
def lp_feas(R, kap, extra_ub=None):
    Aub = [-R, -np.array(LAM3, float)[None, :]]
    r = linprog(np.zeros(8), A_ub=np.vstack(Aub), b_ub=np.zeros(len(R) + 1), A_eq=np.array([kap]), b_eq=[1.0],
                bounds=[(None, None)] * 8, method='highs')
    return r
