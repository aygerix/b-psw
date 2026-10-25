#!/usr/bin/env python3
"""src/cs/e3b/advchk/pen_inspect.py -- inspect Sum m >= 6 pencil points from sq_enum --pencil output: coefficients,
Lucas points per member, symmetry (r_4 = r_5 = -r_2?), shared nodes, reality.  Recomputes everything from (c, l, r).
Usage: pen_inspect.py file [minsum]"""
import sys, numpy as np
fn = sys.argv[1]; minsum = int(sys.argv[2]) if len(sys.argv) > 2 else 6

def fam(z, pencil=(2, 3)):
    c = z[0]; l = np.array(z[1:5]); rr = list(z[5:8])
    r = [rr[0], rr[0], rr[1], rr[2]]          # r_3 := r_2
    return c, l, np.array(r)

def f(c, l, r, j, x):
    A = x**2 - 1
    return 1 + A * (c if j == 1 else l[j - 2] * (x - r[j - 2]))

def table(c, l, r):
    out = []
    for j in range(1, 6):
        coef = [c, 0, 2 - c] if j == 1 else [l[j - 2], -l[j - 2] * r[j - 2], -l[j - 2], l[j - 2] * r[j - 2] + 2]
        pts = []
        for x in np.roots(coef):
            vals = [f(c, l, r, i, x) for i in range(1, 6)]
            if abs(np.prod(vals) + 1) / (np.prod(np.abs(vals)) + 1) < 1e-8:
                pts.append(x)
        out.append(pts)
    return out

seen = []
for line in open(fn):
    st, gd, dg, real, z = eval(line)
    if sum(gd) < minsum:
        continue
    c, l, r = fam(np.array(z))
    key = np.round(np.sort_complex(np.concatenate([[c], l * r, r**2])), 5)
    if any(np.allclose(key, k, atol=1e-5) for k in seen):
        continue
    seen.append(key)
    tab = table(c, l, r)
    sym = abs(r[2] + r[0]) < 1e-6 and abs(r[3] + r[0]) < 1e-6
    sym2 = abs(r[2] + r[3]) < 1e-6
    print(f"[{st}] gd {[len(p) for p in tab]} real={real} c={np.round(c, 6)}")
    print(f"   l = {np.round(l, 6).tolist()}\n   r = {np.round(r, 6).tolist()}  r4=r5=-r2: {sym}  r4=-r5: {sym2}")
    print(f"   Lucas pts: {[np.round(p, 6).tolist() for p in tab]}")
print("distinct Sum >=", minsum, "pencil families:", len(seen))
