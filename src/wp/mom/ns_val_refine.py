import sys, numpy as np
sys.path.insert(0, '.')
import ns_refine as nr
import mpmath as mp
b = np.array([0.517715, 1.990854, 21.097857, 0.896437, 6.526845]); d = np.array([-1.705746, -3.59832, -6.677547, -2.316047, -6.113322])
A, C, hist, smin = nr.refine(list(b), list(b + d), 0, '-', True, 50, 30, verbose=False)
print("half '-' refine: residual history", [mp.nstr(h, 3) for h in hist[:3]], "...", mp.nstr(hist[-1], 3), " smin", mp.nstr(smin, 5))
print("a =", [mp.nstr(z.real, 15) for z in A]); print("c =", [mp.nstr(z.real, 15) for z in C])
r = nr.residual(A, C, 0, '+')
print("'+' residual of the refined half solution:", mp.nstr(max(abs(z) for z in r), 5))
