import numpy as np
from pc7_cert_lib import load, at_shape
D = load()
for R1, R2 in [(-0.9, -0.9), (-0.9, 0.0), (-0.3, 0.3), (0.2, 0.5)]:
    kap, Ac, Bc = at_shape(D, R1, R2)
    print((R1, R2), "kap", np.round(kap, 4))
    print("  nan in A", np.isnan(Ac).sum(), "nan in B", np.isnan(Bc).sum(), " lead", Ac[:, 9], " trail", Ac[:, 1])
