# HD test 1: linear common-root Korselt families f_j = 1 + x/b_j (f_j - 1 = x/b_j, all share root 0 => Korselt automatic).
# Member j Lucas-good iff F(-2 b_j) = -1 iff prod_{i != j} (1 - 2 b_j / b_i) = 1.
# Search real b (distinct, nonzero) with g good members; defect = k - g; HD needs k - g >= k/2.
import numpy as np, itertools, sys
from scipy.optimize import least_squares
rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
k = int(sys.argv[1]); g = int(sys.argv[2])
def res(v, good):
    b = np.concatenate(([1.0], v))
    out = []
    for j in good:
        p = 1.0
        for i in range(k):
            if i != j: p *= (1 - 2*b[j]/b[i])
        out.append(p - 1)
    return np.array(out)
found = []
for good in itertools.combinations(range(k), g):
    for trial in range(400):
        v0 = rng.normal(size=k-1) * rng.choice([0.5, 1, 3])
        try:
            r = least_squares(res, v0, args=(good,), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=2000)
        except Exception: continue
        b = np.concatenate(([1.0], r.x))
        if np.max(np.abs(r.fun)) < 1e-11:
            mind = min(abs(b[i]-b[j]) for i in range(k) for j in range(i))
            if mind > 1e-3 and np.min(np.abs(b)) > 1e-3 and np.max(np.abs(b)) < 1e4:
                found.append((good, b))
    print("good set", good, "solutions so far", len(found), flush=True)
    if found: break
for good, b in found[:10]:
    # count all good members
    allg = [j for j in range(k) if abs(np.prod([1-2*b[j]/b[i] for i in range(k) if i != j]) - 1) < 1e-8]
    print(good, "b =", np.array2string(b, precision=10), "all good:", allg, "signs", np.sign(b))
