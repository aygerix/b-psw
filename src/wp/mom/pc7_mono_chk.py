# check the decomposition against the original evaluator (q22_k7_B.txt) at random data
import numpy as np, sys
exec(open('q22_k7_feas.py').read().split('starts, seed')[0])
from pc7_mono_lib import evalB
rng = np.random.default_rng(1); mx = 0
for _ in range(20):
    r1, r2 = np.sort(rng.uniform(-0.95, 0.95, 2)); w = np.exp(rng.uniform(-2, 1.5, 7))
    a = np.array([eval(c, {}, env(r1, r2, w)) for c in code]); b = evalB(r1, r2, w)
    mx = max(mx, np.max(np.abs(a - b) / (1 + np.abs(a))))
print("max rel diff", mx)
