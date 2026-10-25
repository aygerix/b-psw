"""Boundary systems: m members with w > 0 and 7 - m members at w = 0 (c = infinity), same k = 7 D-equations.
For m = 6 this is 8 equations in 8 unknowns (shape, w_1..w_6).  A nondegenerate positive solution would continue (implicit
function theorem in w_7) to positive solutions of the full 7-member system with small w_7 > 0.
Unknowns r1 = tanh(z0), r2 = tanh(z1), w = exp(z).  Usage: python3 pc7_cert_bnd6.py m starts seed"""
import sys, numpy as np
from scipy.optimize import least_squares
from pc7_cert_dlib import load, Evec
at = load(); m = int(sys.argv[1]); ns = int(sys.argv[2]); rng = np.random.default_rng(int(sys.argv[3]))
def F(z):
    r1, r2 = np.tanh(z[0]), np.tanh(z[1]); ws = np.exp(z[2:]); K, A, C = at(r1, r2)
    return Evec(K, A, C, ws) / Evec(np.abs(K), np.abs(A), np.abs(C), ws)
for s in range(ns):
    z0 = np.concatenate([rng.uniform(-2, 2, 2), rng.uniform(-5, 4, m)])
    try: r = least_squares(F, z0, method='lm' if m + 2 <= 8 else 'trf', max_nfev=3000, xtol=1e-15, ftol=1e-15)
    except Exception as e: print("err", e); continue
    f = np.max(np.abs(F(r.x)))
    tag = "  <-- SOLUTION?" if f < 1e-9 else ""
    print(f"start {s}: scaled res {f:.3e}  r = ({np.tanh(r.x[0]):+.6f},{np.tanh(r.x[1]):+.6f})  w = {np.array2string(np.sort(np.exp(r.x[2:])), precision=5)}{tag}", flush=True)
