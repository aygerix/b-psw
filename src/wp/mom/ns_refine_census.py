#!/usr/bin/env python3
"""ns line: refine every distinct half solution of a census file to high precision (mpmath Gauss-Newton on the
exact '-' remainder system), then report its '+' residual at that precision, the Jacobian conditioning, the reality
margins and the member separation.  Usage: ns_refine_census.py census.npz lin [dps=50] [max=20]"""
import sys, os
import numpy as np
import mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ns_refine as nr

f, lin = sys.argv[1], int(sys.argv[2])
dps = int(sys.argv[3]) if len(sys.argv) > 3 else 50
mx = int(sys.argv[4]) if len(sys.argv) > 4 else 20
d = np.load(f)
order = np.argsort(d['remp'])[:mx]
for j in order:
    a = list(d['a'][j].real); c = list(d['c'][j].real)
    A, C, hist, smin = nr.refine(a, c, lin, '-', True, dps, 40, verbose=False)
    rp = nr.residual(A, C, lin, '+')
    rm = nr.residual(A, C, lin, '-')
    Ar = [z.real for z in A]; Cr = [z.real for z in C]
    marg = min((mp.sqrt(8 * x) - max(abs(y), abs(x - y))) / mp.sqrt(8 * x) for x, y in zip(Ar, Cr))
    sep = min(abs(Ar[i] - Ar[l]) + abs(Cr[i] - Cr[l]) for i in range(len(Ar)) for l in range(i))
    print(f"half solution #{j}: '-' residual {mp.nstr(max(abs(z) for z in rm), 3)} after {len(hist)} its; Jacobian smin "
          f"{mp.nstr(smin, 4)}; '+' residual {mp.nstr(max(abs(z) for z in rp), 6)}; margin {mp.nstr(marg, 6)}; sep {mp.nstr(sep, 6)}")
    print("   a =", [mp.nstr(z, 25) for z in Ar]); print("   c =", [mp.nstr(z, 25) for z in Cr])
