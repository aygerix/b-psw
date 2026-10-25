"""Sign patterns of rho-(c), rho+(c) on the reality interval I_c (see pc7_mono_ext.py), for sample shapes."""
import numpy as np, sys
from pc7_mono_ext_lib import rho
rng = np.random.default_rng(1)
pats = {}
for it in range(int(sys.argv[1])):
    r1, r2 = np.sort(rng.uniform(-1, 1, 2))
    cs = np.exp(np.linspace(-6, 6, 1200)); C = [(c,) + rho(c, r1, r2) for c in cs]
    C = [(c, a, b) for c, ok, a, b in C if ok]
    if not C: continue
    sm = ''.join('+' if a > 0 else '-' for c, a, b in C); sp = ''.join('+' if b > 0 else '-' for c, a, b in C)
    comp = lambda s: ''.join(ch for i, ch in enumerate(s) if i == 0 or s[i - 1] != ch)
    key = (comp(sm), comp(sp)); pats[key] = pats.get(key, 0) + 1
    if pats[key] <= 2: print(f"r = ({r1:+.3f},{r2:+.3f}) e1 = {r1 + r2:+.3f}  I_c = [{C[0][0]:.4g}, {C[-1][0]:.4g}]  rho- {comp(sm)}  rho+ {comp(sp)}")
print(pats)
