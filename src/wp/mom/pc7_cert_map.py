"""Map of the pair relaxation over the shape triangle: max-margin alpha for NORMALISED functionals
  E1/g1^6, E2/g2^6, E3/g1^5, E4/g2^5, E5/h+^6, E6/h-^6, E7/h+^5, E8/h-^5   (g_i = r_i^2 - 1, h+- = H(+-1)).
maximise t: pi_alpha(w, v) >= t * N(w, v) on the (w, v) grid, kap_alpha = 1, lam3_alpha >= 0, |alpha| <= 50.
Usage: python3 pc7_cert_map.py out.npz  (shape list built inside)"""
import sys, time, numpy as np
from scipy.optimize import linprog
from pc7_cert_lib import load, at_shape, LAM3
D = load()
PW = np.exp(np.linspace(np.log(1e-6), np.log(1e6), 49))
WV = np.array([(a, b) for a in PW for b in PW if b >= a])
P = np.array([a ** np.arange(10) for a in PW])                      # (49, 10)
iu = [(i, j) for i in range(len(PW)) for j in range(len(PW)) if j >= i]
def scales(R1, R2):
    g1, g2 = R1**2 - 1, R2**2 - 1; hp, hm = (1 - R1) * (1 - R2), (1 + R1) * (1 + R2)
    return np.array([g1**6, g2**6, g1**5, g2**5, hp**6, hm**6, hp**5, hm**5])
def pairrows(R1, R2):
    kap, Ac, Bc = at_shape(D, R1, R2); s = scales(R1, R2)
    kap = kap / s; Ac = Ac / s[:, None]; Bc = Bc / s[:, None, None]
    Aw = P @ Ac.T                                                     # (49, 8) a_i(w)
    Bw = np.einsum('xm,imn,yn->xyi', P, Bc, P)                        # (49, 49, 8)
    R = np.array([kap / 21 + (Aw[i] + Aw[j]) / 6 + Bw[i, j] for i, j in iu])
    return kap, R
def maxmargin(R1, R2):
    kap, R = pairrows(R1, R2); N = np.abs(R).sum(axis=1, keepdims=True); l3 = LAM3 / scales(R1, R2)
    # vars alpha+ (8), alpha- (8), t ; maximise t ; sum(alpha+ + alpha-) <= 1 ; pi >= t N ; lam3 >= 0 ; kap >= 0
    Rn = R / N
    Aub = np.vstack([np.hstack([-Rn, Rn, np.ones((len(R), 1))]),
                     np.concatenate([-l3, l3, [0]])[None, :], np.concatenate([-kap, kap, [0]])[None, :],
                     np.concatenate([np.ones(16), [0]])[None, :]])
    bub = np.zeros(len(Aub)); bub[-1] = 1.0
    r = linprog(np.concatenate([np.zeros(16), [-1.0]]), A_ub=Aub, b_ub=bub, bounds=[(0, None)] * 16 + [(None, 1)], method='highs')
    if r.status != 0: return (np.nan, np.full(8, np.nan))
    a = r.x[:8] - r.x[8:16]; return (r.x[16], a / np.abs(a).max())
vals = [-0.99, -0.97, -0.93, -0.85, -0.7, -0.5, -0.3, -0.1, 0.0, 0.1, 0.3, 0.5, 0.7, 0.85, 0.93, 0.97, 0.99]
shapes = [(a, b) for a in vals for b in vals if b >= a and a + b >= -1e-12]   # e1 >= 0 suffices by t -> -t
t0 = time.time(); out = []
for s in shapes:
    m, a = maxmargin(*s); out.append((s[0], s[1], m, *a))
    print(f"r = ({s[0]:+.2f},{s[1]:+.2f})  margin {m:+.3e}  alpha {np.round(a, 3)}", flush=True)
np.savez(sys.argv[1], out=np.array(out)); print(f"done {len(shapes)} shapes in {time.time() - t0:.0f}s")
