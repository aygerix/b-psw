#!/usr/bin/env python3
"""src/cs/ident/ctl_check.py -- cross-check the GP controls with the exact Python checker src/bpsw.py.
Each control: (label, factors). Expected: Fermat-side controls pass strong_fermat2 only, Lucas-side controls pass
strong_lucas only; nothing passes both."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from bpsw import strong_fermat2, strong_lucas, bpsw_pass, selfridge_D
from sympy import isprime, fibonacci
ctl = [
 ("Chernick F-form S={1,2,3}, m=6*276", [1+d*6*276 for d in (1,2,3)]),
 ("Chernick F-form S={1,2,3}, m=6*370", [1+d*6*370 for d in (1,2,3)]),
 ("Chernick F-form S={1,2,3,6}, m=6*506", [1+d*6*506 for d in (1,2,3,6)]),
 ("Chernick F-form S={1,4,10,15,30}, m=60*947", [1+d*60*947 for d in (1,4,10,15,30)]),
 ("Lucas-Carmichael L-form S={1,2,3}, m=6*14", [d*6*14-1 for d in (1,2,3)]),
 ("Lucas-Carmichael L-form S={1,2,3}, m=6*70", [d*6*70-1 for d in (1,2,3)]),
 ("Lucas-Carmichael L-form S={1,2,3}, m=6*490", [d*6*490-1 for d in (1,2,3)]),
]
for q in (37, 53, 67):
    F = int(fibonacci(q)); ctl.append((f"composite Fibonacci F_{q}", [F]))
for q in (11, 23, 29):
    ctl.append((f"composite Mersenne M_{q}", [2**q-1]))
for lab, fs in ctl:
    n = 1
    for f in fs: n *= f
    pr = all(isprime(f) for f in fs) if len(fs) > 1 else (not isprime(n))
    print(f"{lab}: n={n}  factors prime/composite-ok={pr}  D*={selfridge_D(n)}  strong_fermat2={strong_fermat2(n)}  strong_lucas={strong_lucas(n)}  BPSW={bpsw_pass(n)}")
