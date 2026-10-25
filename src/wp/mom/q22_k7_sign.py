"""Sign survey of the k = 7 pencil boundary expressions B1..B8 (and B_i +- B_j) over random valid data:
r1, r2 in (-1, 1), c_j > 0 (w_j = 1/c_j > 0). A combination with constant strict sign on the whole domain proves infeasibility."""
import sys, numpy as np
exec(open('q22_k7_feas.py').read().split('starts, seed')[0])   # reuse parser, env()
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 3)
M = int(sys.argv[1])
vals = []
for _ in range(M):
    r1, r2 = rng.uniform(-1, 1, 2)
    w = np.exp(rng.uniform(-6, 4, k))
    d = env(r1, r2, w)
    vals.append([eval(c, {}, d) for c in code])
V = np.array(vals)
names = [f"B{i+1}" for i in range(8)]
combos = {n: V[:, i] for i, n in enumerate(names)}
for i in range(8):
    for j in range(i + 1, 8):
        combos[f"B{i+1}+B{j+1}"] = V[:, i] + V[:, j]; combos[f"B{i+1}-B{j+1}"] = V[:, i] - V[:, j]
for n, v in combos.items():
    pos, neg = np.mean(v > 0), np.mean(v < 0)
    if pos == 1.0 or neg == 1.0: print(f"CONSTANT SIGN {n}: {'+' if pos == 1 else '-'}  min|.| {np.min(np.abs(v)):.3e}")
print("done", M)
