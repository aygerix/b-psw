# HD test 1b: linear common-root families with ALL slopes positive: f_i = 1 + a_i x, a_i = exp(y_i) > 0.
# Good_j: prod_{i != j}(1 - 2 a_i / a_j) = 1.  Look for g good members (HD needs g <= floor(k/2)).
# Also mixed-sign search with sign pattern recorded.
import numpy as np, itertools, sys
from scipy.optimize import least_squares
k = int(sys.argv[1]); g = int(sys.argv[2]); mode = sys.argv[3]; ntr = int(sys.argv[4])
rng = np.random.default_rng(7)
def avec(v):
    return np.concatenate(([1.0], np.exp(v))) if mode == 'pos' else np.concatenate(([1.0], v))
def res(v, good):
    a = avec(v)
    return np.array([np.prod([1 - 2*a[i]/a[j] for i in range(k) if i != j]) - 1 for j in good])
best = {}
nsol = 0; pats = {}
for trial in range(ntr):
    good = tuple(sorted(rng.choice(k, g, replace=False)))
    v0 = rng.normal(size=k-1) * (2 if mode == 'pos' else 1.5)
    try: r = least_squares(res, v0, args=(good,), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
    except Exception: continue
    a = avec(r.x)
    if not np.all(np.isfinite(a)): continue
    m = np.max(np.abs(r.fun))
    md = min(abs(a[i]-a[j]) / max(abs(a[i]), abs(a[j])) for i in range(k) for j in range(i))
    key = good
    if key not in best or m < best[key][0]: best[key] = (m, md, a)
    if m < 1e-12 and md > 1e-4 and np.max(np.abs(a))/np.min(np.abs(a)) < 1e6:
        nsol += 1
        order = np.argsort(a); rank_good = tuple(sorted(int(np.where(order == j)[0][0]) for j in good))
        sp = ''.join('+' if x > 0 else '-' for x in a[order])
        pats[(sp, rank_good)] = pats.get((sp, rank_good), 0) + 1
        if nsol <= 3: print("SOL a =", np.array2string(a, precision=12), "good", good, "res", m, flush=True)
print("k", k, "g", g, mode, "solutions", nsol)
for p, c in sorted(pats.items(), key=lambda t: -t[1])[:20]: print("  sorted sign pattern / good ranks", p, c)
if nsol == 0:
    for key, (m, md, a) in sorted(best.items(), key=lambda t: t[1][0])[:5]:
        print("  best residual for good", key, m, "min rel dist", md, "a", np.array2string(np.sort(a), precision=6))
