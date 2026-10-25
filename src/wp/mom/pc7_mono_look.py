import numpy as np
from pc7_mono_lib import parts
np.set_printoptions(linewidth=200, precision=5)
for r1, r2 in [(-0.3, 0.6), (0.0, 0.2), (-0.5, 0.5)]:
    C0, P1, A, kap = map(np.array, parts(r1, r2))
    print(f"shape ({r1},{r2}) e1={r1+r2:.3f}")
    for i in range(8):
        print(f" B{i+1}: C0 {C0[i]:+.4e} kap {kap[i]:+.0f} phi1 {P1[i][1:]}  A {A[i]}")
