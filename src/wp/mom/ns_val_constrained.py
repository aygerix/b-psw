import sys, numpy as np
sys.path.insert(0, '.')
import ns_core3 as n3
# notebook real '-' solution (6 digits) and a complex one (start 0 of hv5_cx.out)
b = np.array([0.517715, 1.990854, 21.097857, 0.896437, 6.526845]); d = np.array([-1.705746, -3.59832, -6.677547, -2.316047, -6.113322])
bc = np.array([(1.184022-0.657774j), (0.596619+0.071759j), (2.6781+1.020047j), (2.719044-0.722696j), (1.244752+0.201818j)])
dc = np.array([(-3.938373+0.380873j), (-1.976173-0.099753j), (-5.100612-4.548339j), (-5.068403+4.483133j), (-3.199651-0.210308j)])
y = n3.eval_points(5)
for nm, B_, D_, real in (("real", b, d, True), ("complex", bc, dc, False)):
    a = (B_ + 0j)[None]; c = (B_ + D_ + 0j)[None]
    print(nm, "'-' constrained:", np.abs(n3.resid_c(a, c, 0, y, '-', real)).max(), " free:", np.abs(n3.resid(a, c, 0, y, '-', real)).max(),
          " '+' constrained:", np.abs(n3.resid_c(a, c, 0, y, '+', real)).max())
# reflection consistency of the constrained '+' side: '+' residual of x vs '-' residual of reflected x (zero sets agree);
# check on a k = 9 random point that both are O(1) and on a planted-like exact structure we cannot; instead check
# that the fixed part satisfies V(-1), V'(-1):
k = 9
for lin in (0, 1):
    T = np.array([-1 + 1e-6, -1 - 1e-6, -1.0])
    v = n3._vfix(k, lin, T)
    print("lin", lin, "V(-1) =", v[2], "target", 2.0 ** (1 - k), " V'(-1) ~", (v[0] - v[1]) / 2e-6, "target", (1 - k) * 2.0 ** (1 - k))
print("reflection check:")
for nm, B_, D_, real in (("real", b, d, True), ("complex", bc, dc, False)):
    a = (B_ + 0j)[None]; c = (B_ + D_ + 0j)[None]
    print(nm, "'+' constrained at sigma(x):", np.abs(n3.resid_c(a, a - c, 0, y, '+', real)).max(),
          " '+' free at sigma(x):", np.abs(n3.resid(a, a - c, 0, y, '+', real)).max())
