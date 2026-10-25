# ak_collin_check.py -- numerical sanity check of the collinear theorem (exact proof in the notes):
# random collinear admissible families (k cubic members on a line in the (a,c)-plane); at i = argmin a, either
# Re f_l(w_i) > 1 for all l != i, or Re(-f_l(v_i)) > 1 for all l != i.  Also: ray families + linear member (two-ended argument).
import numpy as np
from ak_core import special, fval
rng = np.random.default_rng(5); bad = 0; tested = 0; badray = 0; testedray = 0
def admissible(a, c): return (a > 0) & (a < 32) & (c**2 < 8*a) & ((a - c)**2 < 8*a)
for trial in range(int(__import__("sys").argv[1])):
    k = rng.integers(2, 12)
    # random line through an admissible point
    a0 = np.exp(rng.uniform(np.log(0.01), np.log(31.9))); c0 = a0/2 + (np.sqrt(8*a0) - a0/2)*rng.uniform(-1, 1)
    th = rng.uniform(0, np.pi); ts = rng.normal(0, 1, k)*rng.choice([0.01, 0.3, 3.0])
    a = a0 + ts*np.cos(th); c = c0 + ts*np.sin(th)
    m = admissible(a, c)
    if m.sum() < 2: continue
    a, c = a[m], c[m]; tested += 1
    i = np.argmin(a); w, v = special(a[i:i+1], c[i:i+1]); w, v = w[0], v[0]
    oth = np.arange(len(a)) != i
    okw = np.all(fval(a[oth], c[oth], w).real > 1); okv = np.all((-fval(a[oth], c[oth], v)).real > 1)
    if not (okw or okv): bad += 1
    # ray family (c = u a) + linear member
    u = rng.uniform(-3, 3); aa = np.exp(rng.uniform(np.log(0.01), np.log(31.9), k)); cc = u*aa
    mm = admissible(aa, cc)
    if mm.sum() >= 2:
        aa, cc = aa[mm], cc[mm]; testedray += 1
        found = False
        for i in range(len(aa)):
            w, v = special(aa[i:i+1], cc[i:i+1]); w, v = w[0], v[0]; oth = np.arange(len(aa)) != i
            for z in (w, v):
                vals = np.concatenate([np.abs(fval(aa[oth], cc[oth], z)), [abs(2*z + 1)]])
                pr = np.prod(vals)
                if np.all(vals >= 1) and pr > 1: found = True
                if np.all(vals <= 1) and pr < 1: found = True
        if not found: badray += 1
print("collinear cubic families tested: %d, failures of the argmin-a criterion: %d" % (tested, bad))
print("ray + linear families tested: %d, with no all->=1 or all-<=1 special point: %d" % (testedray, badray))
