"""Write the obstructing configurations (dual of the max-margin LP) at a shape as rationals for an exact gp check.
Usage: python3 pc7_cert_farkas.py nconf r1num/r1den r2num/r2den out.gp"""
import sys, numpy as np
from fractions import Fraction as Fr
from scipy.optimize import linprog
from pc7_cert_dlib import load, Evec
at = load(); rng = np.random.default_rng(5); nconf = int(sys.argv[1]); R1, R2 = Fr(sys.argv[2]), Fr(sys.argv[3])
CF = [np.exp(rng.uniform(-14, 14, 7)) for _ in range(nconf)]
for _ in range(nconf // 2):
    m = rng.integers(1, 4); vals = np.exp(rng.uniform(-14, 14, m)); idx = rng.integers(0, m, 7)
    CF.append(vals[idx] * np.exp(rng.normal(0, 1e-3, 7)))
K, A, C = at(float(R1), float(R2))
R = np.array([Evec(K, A, C, c) for c in CF]); Rn = R / np.abs(R).sum(axis=1, keepdims=True)
r = linprog(np.concatenate([np.zeros(8), [-1]]), A_ub=np.hstack([-Rn, np.ones((len(Rn), 1))]), b_ub=np.zeros(len(Rn)), bounds=[(-1, 1)] * 8 + [(None, None)], method='highs')
du = -r.ineqlin.marginals; order = [i for i in np.argsort(-du) if du[i] > 1e-14]
print("margin", r.x[-1], " #positive duals", len(order))
def rat(x):
    e = int(np.floor(np.log10(x))); return Fr(int(round(x / 10.0**(e - 5))), 1) * Fr(10)**(e - 5)
with open(sys.argv[4], 'w') as f:
    f.write(f"SHAPE = [{R1}, {R2}];\nCONF = [")
    f.write(";".join(",".join(str(rat(x)) for x in CF[i]) for i in order))
    f.write("];\nTHETA = [" + ",".join(f"{du[i]:.6e}" for i in order) + "];\n")
