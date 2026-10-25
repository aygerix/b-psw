"""Advisor's lower bound (13) for E[2^Delta] in W4: prod over ALL levels l = (q, e) of 1 + (2^c_l - 1) 2^-|R_l|,
c_l = bits(q, e) - bits(q, e - 1), R_l = primes of S* needing level >= e at q. Also the exact mean E[Delta]."""
import sys, math
sys.path.insert(0, '.')
from profile2 import primes_upto, load, bits
P = primes_upto(2000)
rows = load(sys.argv[1], P)
ms = [m for _, m in rows]
cnt = {}
for m in ms:
    for q, e in m.items():
        for ee in range(1, e + 1): cnt[(q, ee)] = cnt.get((q, ee), 0) + 1
lb = 0.0; mean = 0.0; top_only = 0.0
em = {}
for m in ms:
    for q, e in m.items(): em[q] = max(em.get(q, 0), e)
for (q, e), c in cnt.items():
    w = bits(q, e) - bits(q, e - 1)
    lb += math.log2(1 + (2 ** w - 1) * 2.0 ** -c)
    mean += w * 2.0 ** -c
    if e == em[q]: top_only += math.log2(1 + (2 ** w - 1) * 2.0 ** -c)
print(f"N = {len(ms)}, levels = {len(cnt)}; lower bound log2 E[2^Delta] >= {lb:.4f} bits (all levels); "
      f"top levels only {top_only:.4f}; exact E[Delta] = {mean:.10f}")
