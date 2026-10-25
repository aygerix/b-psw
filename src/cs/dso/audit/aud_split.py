"""Split-prime usability facts claimed in E2 (from the audit's recomputed rows):
 - D = 5, X = 1e4: number of split primes, those with gcd(odd ord, odd z) = 1, and whether 5 | ord*z for them;
 - index inequality: for split p with coprime odd orders, odd(i1)*odd(i2) >= odd(p-1) (i1 = (p-1)/ord, i2 = (p-1)/z);
 - usable split primes at X = 2e4, 6e4 (with their orders)."""
import math
v2 = lambda x: (x & -x).bit_length() - 1
odd = lambda x: x >> v2(x)
for X in (10000, 20000, 60000):
    R = [tuple(map(int, l.split())) for l in open(f"rows_X{X}.txt")]
    sp = [r for r in R if r[3] == 1]
    cop = [r for r in sp if math.gcd(odd(r[1]), odd(r[2])) == 1]
    us = [r for r in cop if r[1] % 5 and r[2] % 5]
    viol = [r for r in cop if odd((r[0]-1)//r[1]) * odd((r[0]-1)//r[2]) < odd(r[0]-1)]
    print(f"X = {X}: split {len(sp)}, coprime odd parts {len(cop)} {cop if len(cop) < 8 else ''}, usable (5 not | ord z) {len(us)} {us}; index-inequality violations {len(viol)}")
