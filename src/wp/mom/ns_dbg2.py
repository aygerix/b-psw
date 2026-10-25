import sys, numpy as np
sys.path.insert(0, '.')
import ns_core3 as n3
b = np.array([0.517715, 1.990854, 21.097857, 0.896437, 6.526845])
d = np.array([-1.705746, -3.59832, -6.677547, -2.316047, -6.113322])
a = (b + 0j)[None]; c = (b + d + 0j)[None]
y = n3.eval_points(5)
for j in range(6):
    print("J =", j, " partial residual:", np.abs(n3.resid_partial(a, c, 0, y, 5, j, True)).max())
print("full '-' (fixed-lead W):", np.abs(n3.resid(a, c, 0, y, '-', True)).max())
