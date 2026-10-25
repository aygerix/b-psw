"""Sign survey of the evaluation-form boundary expressions:
 EV1 = V(r1) - c0 G(r1)^6 = B1 + B2 r1, EV2 = B1 + B2 r2, DV1 = V'(r1) - c0 (G^6)'(r1) = B3 + B4 r1, DV2,
 EW+ = W(1) - c0 H(1)^6 = B5 + B6, EW- = B5 - B6, DW+ = B7 + B8, DW- = B7 - B8,
over random shapes and 7-member w, both on (0, inf) and on the reality interval I."""
import numpy as np, sys
from pc7_mono_lib import parts
from pc7_mono_lpI import interval
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
NS = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
def fastB(pr, w):
    C0, P1, A, kap = pr
    pw = np.array([np.sum(w**a) for a in range(10)])
    e3 = (pw[1]**3 - 3 * pw[1] * pw[2] + 2 * pw[3]) / 6
    # pair sums S_n = sum_{j<l} u v h_n(u, v)
    S = np.zeros(5)
    for j in range(len(w)):
        u = w[j]; v = w[j + 1:]
        for n in range(5): S[n] += np.sum(u * v * sum(u**m * v**(n - m) for m in range(n + 1)))
    return np.array([C0[i] + P1[i][1:] @ pw[1:] + A[i] @ S + kap[i] * e3 for i in range(8)])
names = ['EV1', 'EV2', 'DV1', 'DV2', 'EW+', 'EW-', 'DW+', 'DW-']
cnt = {m: {nm: [0, 0] for nm in names} for m in ('pos', 'I')}
for _ in range(NS):
    r1, r2 = np.sort(rng.uniform(-1, 1, 2))
    I = interval(r1, r2)
    pr = tuple(map(np.array, parts(r1, r2)))
    for mode in ('pos', 'I'):
        if mode == 'pos': w = np.exp(rng.uniform(-6, 4, 7))
        else:
            if I is None: continue
            lo, hi = I; hi = min(hi, 1e4 * lo); w = np.exp(rng.uniform(np.log(lo), np.log(hi), 7))
        B = fastB(pr, w)
        E = [B[0] + B[1] * r1, B[0] + B[1] * r2, B[2] + B[3] * r1, B[2] + B[3] * r2, B[4] + B[5], B[4] - B[5], B[6] + B[7], B[6] - B[7]]
        for nm, e in zip(names, E): cnt[mode][nm][0 if e > 0 else 1] += 1
for mode in cnt: print(mode, {nm: tuple(v) for nm, v in cnt[mode].items()})
