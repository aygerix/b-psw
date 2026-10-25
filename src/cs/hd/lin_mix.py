# HD test 1d: mixed-sign linear common-root k = 5: can 4 members be Lucas-good (defect 1)? 4 equations, 4 unknowns (a1 = 1).
# Also k = 3 mixed, 2 good (shown impossible by hand), and k = 7 mixed, 6 good.
import numpy as np, sys, time, itertools
from scipy.optimize import least_squares
k = int(sys.argv[1]); g = int(sys.argv[2]); ntr = int(sys.argv[3]); seed = int(sys.argv[4])
rng = np.random.default_rng(seed)
def res(v, good):
    a = np.concatenate(([1.0], v))
    return np.array([np.prod([1 - 2*a[i]/a[j] for i in range(k) if i != j]) - 1 for j in good])
t0 = time.time(); sols = []
for tr in range(ntr):
    good = sorted(rng.choice(k, g, replace=False))
    v0 = rng.normal(size=k-1) * rng.choice([0.5, 1, 3])
    try: r = least_squares(res, v0, args=(good,), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
    except Exception: continue
    a = np.concatenate(([1.0], r.x))
    if not np.all(np.isfinite(a)): continue
    md = min(abs(a[i]-a[j])/max(abs(a[i]), abs(a[j])) for i in range(k) for j in range(i))
    if np.max(np.abs(r.fun)) < 1e-12 and md > 1e-6 and np.max(np.abs(a))/np.min(np.abs(a)) < 1e7:
        sols.append(a)
        if len(sols) <= 5: print("SOL", np.array2string(np.sort(a), precision=12, max_line_width=200), "good(unsorted idx)", good, flush=True)
print(f"k={k} g={g} trials={ntr} solutions={len(sols)} time={time.time()-t0:.1f}s")
