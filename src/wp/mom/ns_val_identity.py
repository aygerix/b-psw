import sys, numpy as np
sys.path.insert(0, '.')
import ns_core3 as n3
b = np.array([0.517715, 1.990854, 21.097857, 0.896437, 6.526845])
d = np.array([-1.705746, -3.59832, -6.677547, -2.316047, -6.113322])
a = (b + 0j)[None]; c = (b + d + 0j)[None]
y = n3.eval_points(5)
print("notebook '-' solution, identity residual '-':", np.abs(n3.resid(a, c, 0, y, '-', True)).max(),
      " '+':", np.abs(n3.resid(a, c, 0, y, '+', True)).max())
rs = np.random.default_rng(1)
a = rs.uniform(0.5, 3, (3, 5)) + 0j; c = rs.uniform(-1, 1, (3, 5)) + 0j
print("random point, '-':", np.abs(n3.resid(a, c, 0, y, '-', True)).max(axis=1))
# diagonal: all members equal -> NOT a solution
a = np.full((1, 5), 1.3) + 0j; c = np.full((1, 5), 0.4) + 0j
print("diagonal, '-+':", np.abs(n3.resid(a, c, 0, y, '-+', True)).max())
