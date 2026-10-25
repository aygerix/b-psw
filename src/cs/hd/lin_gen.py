# HD test 1e: positive linear common-root families, ordered a_1 = 1 < a_2 < ... < a_k; impose g good members chosen at
# random among the interior ranks 2..k-1 (ranks 1, k are always dead: sign lemma). Log every solution found.
import numpy as np, sys, time
from scipy.optimize import least_squares
k = int(sys.argv[1]); g = int(sys.argv[2]); ntr = int(sys.argv[3]); seed = int(sys.argv[4])
rng = np.random.default_rng(seed)
def avec(y):
    a = [1.0]
    for t in y: a.append(a[-1] * (1 + np.exp(t)))
    return np.array(a)
def res(y, good):
    a = avec(y)
    return np.array([np.prod([1 - 2*a[i]/a[j] for i in range(k) if i != j]) - 1 for j in good])
t0 = time.time(); ns = 0; pats = {}
for tr in range(ntr):
    good = sorted(rng.choice(np.arange(1, k-1), g, replace=False))
    y0 = rng.normal(size=k-1) * rng.choice([0.5, 1.0, 2.0]) + rng.choice([-1.5, -0.5, 0.0, 1.0])
    try: r = least_squares(res, y0, args=(good,), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
    except Exception: continue
    a = avec(r.x)
    if np.all(np.isfinite(a)) and np.max(np.abs(r.fun)) < 1e-12 and np.min(np.diff(a)/a[1:]) > 1e-5 and a[-1] < 1e8:
        ns += 1; pats[tuple(good)] = pats.get(tuple(good), 0) + 1
        print("SOL good_ranks", [int(x) for x in good], "a", np.array2string(a, precision=13, max_line_width=250), flush=True)
print(f"k={k} g={g} trials={ntr} solutions={ns} patterns={pats} time={time.time()-t0:.1f}s", flush=True)
