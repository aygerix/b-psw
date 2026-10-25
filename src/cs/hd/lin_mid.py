# HD test 1c: positive linear common-root families, ordered a_1 < ... < a_k (a_1 = 1, a_{i+1} = a_i (1 + e^{y_i})),
# impose Lucas-good on members in GOOD (rank indices); default: all middle members 2..k-1 (defect exactly 2 = sign-lemma bound).
# Mixed mode: a_i = signed, impose g good members (any).
import numpy as np, sys, time
from scipy.optimize import least_squares
k = int(sys.argv[1]); ntr = int(sys.argv[2]); seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
good = list(range(1, k-1)) if len(sys.argv) <= 4 else [int(t) for t in sys.argv[4].split(',')]
rng = np.random.default_rng(seed)
def avec(y):
    a = [1.0]
    for t in y: a.append(a[-1] * (1 + np.exp(t)))
    return np.array(a)
def res(y):
    a = avec(y)
    # log form for stability: sum log|1 - 2 a_i/a_j| = 0 and sign product = +1
    out = []
    for j in good:
        f = np.array([1 - 2*a[i]/a[j] for i in range(k) if i != j])
        out.append(np.prod(f) - 1)
    return np.array(out)
t0 = time.time(); sols = []
for tr in range(ntr):
    y0 = rng.normal(size=k-1) * rng.choice([0.5, 1.0, 2.0]) + rng.choice([-1.0, 0.0, 1.0])
    try: r = least_squares(res, y0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
    except Exception: continue
    a = avec(r.x)
    if np.all(np.isfinite(a)) and np.max(np.abs(r.fun)) < 1e-12 and np.min(np.diff(a)/a[1:]) > 1e-6 and a[-1] < 1e8:
        sols.append(a)
        if len(sols) <= 4: print("SOL", np.array2string(a, precision=12, max_line_width=200), "maxres", np.max(np.abs(r.fun)), flush=True)
print(f"k={k} good={good} trials={ntr} solutions={len(sols)} time={time.time()-t0:.1f}s")
