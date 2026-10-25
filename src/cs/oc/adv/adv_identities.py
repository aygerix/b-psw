#!/usr/bin/env python3
"""Independent check of the 'largest piece' identities (line OC claim 9), k <= 300, q odd <= 301."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from adv_rows import F, L
bad = 0
for k in range(1, 301):
    t = [(F(4*k+1) - 1, F(2*k) * L(2*k+1)), (F(4*k+3) - 1, F(2*k+2) * L(2*k+1)), (F(4*k+1) + 1, F(2*k+1) * L(2*k)),
         (F(4*k+3) + 1, F(2*k+1) * L(2*k+2)), (L(4*k+1) - 1, 5 * F(2*k) * F(2*k+1)), (L(4*k+3) - 1, L(2*k+1) * L(2*k+2)),
         (L(4*k+1) + 1, L(2*k) * L(2*k+1)), (L(4*k+3) + 1, 5 * F(2*k+1) * F(2*k+2))]
    bad += sum(a != b for a, b in t)
for q in range(5, 302, 2):
    k = (q - 1) // 2
    bad += ((2**q + 1) // 3 + 1) * 3 != 4 * (2**(q-2) + 1)
    bad += 2**q + 2**(k+1) + 1 + 1 != 2 * (4**k + 2**k + 1)
    bad += 2**q - 2**(k+1) + 1 + 1 != 2 * (4**k - 2**k + 1)
    bad += 2**q - 2**(k+1) + 1 - 1 != 2**(k+1) * (2**k - 1)
    bad += 2**q + 2**(k+1) + 1 - 1 != 2**(k+1) * (2**k + 1)
print('identity failures:', bad, '(8 Fibonacci/Lucas identities k<=300; Wagstaff and Gaussian-Mersenne p-+1 forms, odd q<=301)')
