"""Sparsest (L1-minimal) separating weights alpha at a few shapes; also test whether a FIXED small support works everywhere."""
import sys, itertools, numpy as np
from scipy.optimize import linprog
exec(open('q22_k7_feas.py').read().split('starts, seed')[0])
rng = np.random.default_rng(5)
def samples(r1, r2, ns):
    W = np.exp(rng.uniform(-7, 5, (ns, k)))
    Bm = np.array([[eval(c, {}, env(r1, r2, w)) for c in code] for w in W])
    return Bm / np.maximum(1.0, np.abs(Bm).max(axis=1, keepdims=True))
def sep(Bm, supp):
    n = len(supp); A = Bm[:, supp]
    # variables alpha+ , alpha- (2n); minimise sum; constraint A(a+ - a-) >= 1e-6
    c = np.ones(2 * n); Aub = -np.hstack([A, -A]); bub = -1e-6 * np.ones(len(A))
    r = linprog(c, A_ub=Aub, b_ub=bub, bounds=[(0, 1e6)] * (2 * n), method='highs')
    if r.status != 0: return None
    a = r.x[:n] - r.x[n:]; return a / np.max(np.abs(a))
shapes = []
for r1, r2 in shapes:
    Bm = samples(r1, r2, 600); a = sep(Bm, list(range(8)))
    print(f"r = ({r1:+.2f},{r2:+.2f}) alpha =", None if a is None else np.round(a, 4))
# fixed supports of size 2 and 3: which work at ALL shapes?
grid = [(a, b) for a in np.linspace(-0.9, 0.9, 7) for b in np.linspace(-0.9, 0.9, 7) if b > a + 0.05]
data = {s: samples(*s, 300) for s in grid}
for size in (4, 5):
    good = []
    for supp in itertools.combinations(range(8), size):
        if all(sep(data[s], list(supp)) is not None for s in grid): good.append(tuple(i + 1 for i in supp))
    print(f"supports of size {size} separating at all {len(grid)} shapes: {good}")
